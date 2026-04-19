#include "common.hpp"

#include <cmath>
#include <complex>
#include <memory>
#include <utility>
#include <vector>

namespace
{

using mfem::Array;
using mfem::Array2D;
using mfem::BilinearForm;
using mfem::ComplexGridFunction;
using mfem::ComplexLinearForm;
using mfem::ConstantCoefficient;
using mfem::CurlCurlIntegrator;
using mfem::ElementTransformation;
using mfem::Geometry;
using mfem::GridFunction;
using mfem::IntRules;
using mfem::Mesh;
using mfem::RestrictedCoefficient;
using mfem::ScalarVectorProductCoefficient;
using mfem::SesquilinearForm;
using mfem::Vector;
using mfem::VectorCoefficient;
using mfem::VectorFEDomainLFIntegrator;
using mfem::VectorFEMassIntegrator;
using mfem::VectorFunctionCoefficient;
using mfem::VectorRestrictedCoefficient;

enum ProbType
{
   beam = 0,
   disc = 1,
   lshape = 2,
   fichera = 3,
   load_src = 4
};

template <typename T> T pow2(const T &x) { return x * x; }

class PmlRegion
{
private:
   Mesh *mesh;
   int dim;
   Array2D<mfem::real_t> length;
   Array2D<mfem::real_t> comp_dom_bdr;
   Array2D<mfem::real_t> dom_bdr;
   Array<int> elems;

   void SetBoundaries()
   {
      comp_dom_bdr.SetSize(dim, 2);
      dom_bdr.SetSize(dim, 2);
      Vector pmin;
      Vector pmax;
      mesh->GetBoundingBox(pmin, pmax);
      for (int i = 0; i < dim; i++)
      {
         dom_bdr(i, 0) = pmin(i);
         dom_bdr(i, 1) = pmax(i);
         comp_dom_bdr(i, 0) = dom_bdr(i, 0) + length(i, 0);
         comp_dom_bdr(i, 1) = dom_bdr(i, 1) - length(i, 1);
      }
   }

public:
   PmlRegion(Mesh *mesh_, const Array2D<mfem::real_t> &length_)
      : mesh(mesh_), dim(mesh_->Dimension()), length(length_)
   {
      SetBoundaries();
   }

   const Array2D<mfem::real_t> &GetCompDomainBdr() const { return comp_dom_bdr; }
   const Array<int> &GetMarkedPMLElements() const { return elems; }

   void SetAttributes(Mesh *mesh_)
   {
      for (int i = 0; i < mesh_->GetNBE(); ++i)
      {
         mesh_->GetBdrElement(i)->SetAttribute(i + 1);
      }

      const int nrelem = mesh_->GetNE();
      elems.SetSize(nrelem);

      for (int i = 0; i < nrelem; ++i)
      {
         elems[i] = 1;
         bool in_pml = false;
         mfem::Element *el = mesh_->GetElement(i);
         Array<int> vertices;
         el->SetAttribute(1);
         el->GetVertices(vertices);

         for (int iv = 0; iv < vertices.Size(); ++iv)
         {
            mfem::real_t *coords = mesh_->GetVertex(vertices[iv]);
            for (int comp = 0; comp < dim; ++comp)
            {
               if (coords[comp] > comp_dom_bdr(comp, 1) ||
                   coords[comp] < comp_dom_bdr(comp, 0))
               {
                  in_pml = true;
                  break;
               }
            }
            if (in_pml) { break; }
         }

         if (in_pml)
         {
            elems[i] = 0;
            el->SetAttribute(2);
         }
      }
      mesh_->SetAttributes();
   }

   void StretchFunction(const Vector &x,
                        std::vector<std::complex<mfem::real_t>> &dxs,
                        mfem::real_t mu,
                        mfem::real_t epsilon,
                        mfem::real_t omega) const
   {
      constexpr std::complex<mfem::real_t> zi(0.0, 1.0);
      const mfem::real_t n = 2.0;
      const mfem::real_t c = 5.0;
      const mfem::real_t k = omega * std::sqrt(epsilon * mu);

      for (int i = 0; i < dim; ++i)
      {
         dxs[i] = 1.0;
         if (x(i) >= comp_dom_bdr(i, 1))
         {
            const mfem::real_t coeff = n * c / k / std::pow(length(i, 1), n);
            dxs[i] = 1.0 + zi * coeff *
                     std::abs(std::pow(x(i) - comp_dom_bdr(i, 1), n - 1.0));
         }
         if (x(i) <= comp_dom_bdr(i, 0))
         {
            const mfem::real_t coeff = n * c / k / std::pow(length(i, 0), n);
            dxs[i] = 1.0 + zi * coeff *
                     std::abs(std::pow(x(i) - comp_dom_bdr(i, 0), n - 1.0));
         }
      }
   }
};

class Ex25PmlContext;

class PmlDiagMatrixCoefficient : public VectorCoefficient
{
private:
   Ex25PmlContext *context;
   void (Ex25PmlContext::*function)(const Vector &, Vector &) const;

public:
   PmlDiagMatrixCoefficient(int dim,
                            Ex25PmlContext *ctx,
                            void (Ex25PmlContext::*fn)(const Vector &, Vector &) const)
      : VectorCoefficient(dim), context(ctx), function(fn)
   {}

   using VectorCoefficient::Eval;

   void Eval(Vector &K, ElementTransformation &T,
             const mfem::IntegrationPoint &ip) override;
};

class Ex25PmlContext
{
private:
   Mesh *mesh;
   ProbType prob;
   int dim;
   mfem::real_t mu = 1.0;
   mfem::real_t epsilon = 1.0;
   mfem::real_t omega = 1.0;
   std::unique_ptr<PmlRegion> pml;
   std::vector<std::unique_ptr<mfem::Coefficient>> op_scalar_coeffs;
   std::vector<std::unique_ptr<mfem::VectorCoefficient>> op_vector_coeffs;
   std::vector<std::unique_ptr<mfem::Coefficient>> prec_scalar_coeffs;
   std::vector<std::unique_ptr<mfem::VectorCoefficient>> prec_vector_coeffs;

   template <typename T, typename... Args>
   T &OwnOpScalar(Args &&... args)
   {
      op_scalar_coeffs.push_back(
         std::make_unique<T>(std::forward<Args>(args)...));
      return static_cast<T &>(*op_scalar_coeffs.back());
   }

   template <typename T, typename... Args>
   T &OwnOpVector(Args &&... args)
   {
      op_vector_coeffs.push_back(
         std::make_unique<T>(std::forward<Args>(args)...));
      return static_cast<T &>(*op_vector_coeffs.back());
   }

   template <typename T, typename... Args>
   T &OwnPrecScalar(Args &&... args)
   {
      prec_scalar_coeffs.push_back(
         std::make_unique<T>(std::forward<Args>(args)...));
      return static_cast<T &>(*prec_scalar_coeffs.back());
   }

   template <typename T, typename... Args>
   T &OwnPrecVector(Args &&... args)
   {
      prec_vector_coeffs.push_back(
         std::make_unique<T>(std::forward<Args>(args)...));
      return static_cast<T &>(*prec_vector_coeffs.back());
   }

public:
   Ex25PmlContext(Mesh *mesh_, int prob_int)
      : mesh(mesh_), prob(static_cast<ProbType>(prob_int)), dim(mesh_->Dimension())
   {
      Array2D<mfem::real_t> length(dim, 2);
      length = 0.0;
      switch (prob)
      {
         case disc:
            length = 0.2;
            break;
         case lshape:
            length(0, 0) = 0.1;
            length(1, 0) = 0.1;
            break;
         case fichera:
            length(0, 1) = 0.5;
            length(1, 1) = 0.5;
            length(2, 1) = 0.5;
            break;
         case beam:
            length(0, 1) = 2.0;
            break;
         case load_src:
         default:
            length = 0.25;
            break;
      }
      pml = std::make_unique<PmlRegion>(mesh, length);
   }

   void SetParameters(mfem::real_t mu_, mfem::real_t epsilon_, mfem::real_t omega_)
   {
      mu = mu_;
      epsilon = epsilon_;
      omega = omega_;
   }

   void SetAttributes(Mesh *mesh_) { pml->SetAttributes(mesh_); }

   void BuildEssentialBoundary(const Mesh &mesh_, Array<int> &ess_bdr) const
   {
      ess_bdr = 1;
      if (prob != lshape && prob != fichera) { return; }

      ess_bdr = 0;
      for (int j = 0; j < mesh_.GetNBE(); j++)
      {
         Vector center(dim);
         const int bdrgeom = mesh_.GetBdrElementGeometry(j);
         ElementTransformation *tr =
            const_cast<Mesh &>(mesh_).GetBdrElementTransformation(j);
         tr->Transform(mfem::Geometries.GetCenter(bdrgeom), center);
         const int attr = mesh_.GetBdrAttribute(j);
         switch (prob)
         {
            case lshape:
               if (center[0] == 1.0 || center[0] == 0.5 || center[1] == 0.5)
               {
                  ess_bdr[attr - 1] = 1;
               }
               break;
            case fichera:
               if (center[0] == -1.0 || center[0] == 0.0 ||
                   center[1] == 0.0 || center[2] == 0.0)
               {
                  ess_bdr[attr - 1] = 1;
               }
               break;
            default:
               break;
         }
      }
   }

   void Source(const Vector &x, Vector &f) const
   {
      Vector center(dim);
      mfem::real_t r = 0.0;
      for (int i = 0; i < dim; ++i)
      {
         center(i) = 0.5 * (pml->GetCompDomainBdr()(i, 0) +
                            pml->GetCompDomainBdr()(i, 1));
         r += pow2(x[i] - center[i]);
      }
      const mfem::real_t n = 5.0 * omega * std::sqrt(epsilon * mu) / M_PI;
      const mfem::real_t coeff = pow2(n) / M_PI;
      const mfem::real_t alpha = -pow2(n) * r;
      f = 0.0;
      f[0] = coeff * std::exp(alpha);
   }

   void MaxwellSolution(const Vector &x,
                        std::vector<std::complex<mfem::real_t>> &E) const
   {
      constexpr std::complex<mfem::real_t> zi(0.0, 1.0);
      const mfem::real_t k = omega * std::sqrt(epsilon * mu);

      for (int i = 0; i < dim; ++i) { E[i] = 0.0; }

      switch (prob)
      {
         case disc:
         case lshape:
         case fichera:
         {
            Vector shift(dim);
            shift = 0.0;
            if (prob == fichera) { shift = 1.0; }
            if (prob == disc) { shift = -0.5; }
            if (prob == lshape) { shift = -1.0; }

            if (dim == 2)
            {
               const mfem::real_t x0 = x(0) + shift(0);
               const mfem::real_t x1 = x(1) + shift(1);
               const mfem::real_t r = std::sqrt(x0 * x0 + x1 * x1);
               const mfem::real_t beta = k * r;

               const std::complex<mfem::real_t> Ho =
                  mfem::real_t(jn(0, beta)) + zi * mfem::real_t(yn(0, beta));
               const std::complex<mfem::real_t> Ho_r =
                  -k * (mfem::real_t(jn(1, beta)) + zi * mfem::real_t(yn(1, beta)));
               const std::complex<mfem::real_t> Ho_rr =
                  -k * k * (1.0 / beta *
                            (mfem::real_t(jn(1, beta)) + zi * mfem::real_t(yn(1, beta))) -
                            (mfem::real_t(jn(2, beta)) + zi * mfem::real_t(yn(2, beta))));

               const mfem::real_t r_x = x0 / r;
               const mfem::real_t r_y = x1 / r;
               const mfem::real_t r_xy = -(r_x / r) * r_y;
               const mfem::real_t r_xx = (1.0 / r) * (1.0 - r_x * r_x);

               const std::complex<mfem::real_t> val = 0.25 * zi * Ho;
               const std::complex<mfem::real_t> val_xx =
                  0.25 * zi * (r_xx * Ho_r + r_x * r_x * Ho_rr);
               const std::complex<mfem::real_t> val_xy =
                  0.25 * zi * (r_xy * Ho_r + r_x * r_y * Ho_rr);
               E[0] = zi / k * (k * k * val + val_xx);
               E[1] = zi / k * val_xy;
            }
            else if (dim == 3)
            {
               const mfem::real_t x0 = x(0) + shift(0);
               const mfem::real_t x1 = x(1) + shift(1);
               const mfem::real_t x2 = x(2) + shift(2);
               const mfem::real_t r = std::sqrt(x0 * x0 + x1 * x1 + x2 * x2);

               const mfem::real_t r_x = x0 / r;
               const mfem::real_t r_y = x1 / r;
               const mfem::real_t r_z = x2 / r;
               const mfem::real_t r_xx = (1.0 / r) * (1.0 - r_x * r_x);
               const mfem::real_t r_yx = -(r_y / r) * r_x;
               const mfem::real_t r_zx = -(r_z / r) * r_x;

               const std::complex<mfem::real_t> val = std::exp(zi * k * r) / r;
               const std::complex<mfem::real_t> val_r = val / r * (zi * k * r - 1.0);
               const std::complex<mfem::real_t> val_rr =
                  val / (r * r) * (-k * k * r * r
                                   - 2.0 * zi * k * r + 2.0);

               const std::complex<mfem::real_t> val_xx = val_rr * r_x * r_x +
                                                         val_r * r_xx;
               const std::complex<mfem::real_t> val_yx = val_rr * r_x * r_y +
                                                         val_r * r_yx;
               const std::complex<mfem::real_t> val_zx = val_rr * r_x * r_z +
                                                         val_r * r_zx;
               const std::complex<mfem::real_t> alpha =
                  zi * k / 4.0 / M_PI / k / k;
               E[0] = alpha * (k * k * val + val_xx);
               E[1] = alpha * val_yx;
               E[2] = alpha * val_zx;
            }
            break;
         }
         case beam:
            if (dim == 3)
            {
               const mfem::real_t k10 = std::sqrt(k * k - M_PI * M_PI);
               E[1] = -zi * k / M_PI *
                      std::sin(M_PI * x(2)) * std::exp(zi * k10 * x(0));
            }
            else if (dim == 2)
            {
               E[1] = -zi * k / M_PI * std::exp(zi * k * x(0));
            }
            break;
         case load_src:
         default:
            break;
      }
   }

   void ExactReal(const Vector &x, Vector &E) const
   {
      std::vector<std::complex<mfem::real_t>> eval(E.Size());
      MaxwellSolution(x, eval);
      for (int i = 0; i < dim; ++i) { E[i] = eval[i].real(); }
   }

   void ExactImag(const Vector &x, Vector &E) const
   {
      std::vector<std::complex<mfem::real_t>> eval(E.Size());
      MaxwellSolution(x, eval);
      for (int i = 0; i < dim; ++i) { E[i] = eval[i].imag(); }
   }

   void BoundaryReal(const Vector &x, Vector &E) const
   {
      E = 0.0;
      bool in_pml = false;
      for (int i = 0; i < dim; ++i)
      {
         if (x(i) - pml->GetCompDomainBdr()(i, 0) < 0.0 ||
             x(i) - pml->GetCompDomainBdr()(i, 1) > 0.0)
         {
            in_pml = true;
            break;
         }
      }
      if (!in_pml) { ExactReal(x, E); }
   }

   void BoundaryImag(const Vector &x, Vector &E) const
   {
      E = 0.0;
      bool in_pml = false;
      for (int i = 0; i < dim; ++i)
      {
         if (x(i) - pml->GetCompDomainBdr()(i, 0) < 0.0 ||
             x(i) - pml->GetCompDomainBdr()(i, 1) > 0.0)
         {
            in_pml = true;
            break;
         }
      }
      if (!in_pml) { ExactImag(x, E); }
   }

   void DetJJtJInvRe(const Vector &x, Vector &D) const
   {
      std::vector<std::complex<mfem::real_t>> dxs(dim);
      std::complex<mfem::real_t> det(1.0, 0.0);
      pml->StretchFunction(x, dxs, mu, epsilon, omega);
      for (int i = 0; i < dim; ++i) { det *= dxs[i]; }
      for (int i = 0; i < dim; ++i) { D(i) = (det / pow2(dxs[i])).real(); }
   }

   void DetJJtJInvIm(const Vector &x, Vector &D) const
   {
      std::vector<std::complex<mfem::real_t>> dxs(dim);
      std::complex<mfem::real_t> det(1.0, 0.0);
      pml->StretchFunction(x, dxs, mu, epsilon, omega);
      for (int i = 0; i < dim; ++i) { det *= dxs[i]; }
      for (int i = 0; i < dim; ++i) { D(i) = (det / pow2(dxs[i])).imag(); }
   }

   void DetJJtJInvAbs(const Vector &x, Vector &D) const
   {
      std::vector<std::complex<mfem::real_t>> dxs(dim);
      std::complex<mfem::real_t> det(1.0, 0.0);
      pml->StretchFunction(x, dxs, mu, epsilon, omega);
      for (int i = 0; i < dim; ++i) { det *= dxs[i]; }
      for (int i = 0; i < dim; ++i) { D(i) = std::abs(det / pow2(dxs[i])); }
   }

   void DetJInvJtJRe(const Vector &x, Vector &D) const
   {
      std::vector<std::complex<mfem::real_t>> dxs(dim);
      std::complex<mfem::real_t> det(1.0, 0.0);
      pml->StretchFunction(x, dxs, mu, epsilon, omega);
      for (int i = 0; i < dim; ++i) { det *= dxs[i]; }
      if (dim == 2)
      {
         D = (1.0 / det).real();
      }
      else
      {
         for (int i = 0; i < dim; ++i) { D(i) = (pow2(dxs[i]) / det).real(); }
      }
   }

   void DetJInvJtJIm(const Vector &x, Vector &D) const
   {
      std::vector<std::complex<mfem::real_t>> dxs(dim);
      std::complex<mfem::real_t> det(1.0, 0.0);
      pml->StretchFunction(x, dxs, mu, epsilon, omega);
      for (int i = 0; i < dim; ++i) { det *= dxs[i]; }
      if (dim == 2)
      {
         D = (1.0 / det).imag();
      }
      else
      {
         for (int i = 0; i < dim; ++i) { D(i) = (pow2(dxs[i]) / det).imag(); }
      }
   }

   void DetJInvJtJAbs(const Vector &x, Vector &D) const
   {
      std::vector<std::complex<mfem::real_t>> dxs(dim);
      std::complex<mfem::real_t> det(1.0, 0.0);
      pml->StretchFunction(x, dxs, mu, epsilon, omega);
      for (int i = 0; i < dim; ++i) { det *= dxs[i]; }
      if (dim == 2)
      {
         D = std::abs(1.0 / det);
      }
      else
      {
         for (int i = 0; i < dim; ++i) { D(i) = std::abs(pow2(dxs[i]) / det); }
      }
   }

   void ProjectBoundary(ComplexGridFunction &x, Array<int> &ess_bdr) const
   {
      VectorFunctionCoefficient real_coef(dim,
                                          [this](const Vector &xx, Vector &yy)
      {
         BoundaryReal(xx, yy);
      });
      VectorFunctionCoefficient imag_coef(dim,
                                          [this](const Vector &xx, Vector &yy)
      {
         BoundaryImag(xx, yy);
      });
      x.ProjectBdrCoefficientTangent(real_coef, imag_coef, ess_bdr);
   }

   void AssembleRhs(ComplexLinearForm &b) const
   {
      if (prob == load_src)
      {
         VectorFunctionCoefficient f(dim, [this](const Vector &xx, Vector &yy)
         {
            Source(xx, yy);
         });
         b.AddDomainIntegrator(nullptr, new VectorFEDomainLFIntegrator(f));
      }
      b = 0.0;
      b.Assemble();
   }

   void AddOperatorIntegrators(SesquilinearForm &a)
   {
      Array<int> attr;
      Array<int> attr_pml;
      if (mesh->attributes.Size())
      {
         attr.SetSize(mesh->attributes.Max());
         attr_pml.SetSize(mesh->attributes.Max());
         attr = 0;
         attr[0] = 1;
         attr_pml = 0;
         if (mesh->attributes.Max() > 1) { attr_pml[1] = 1; }
      }

      auto &muinv = OwnOpScalar<ConstantCoefficient>(1.0 / mu);
      auto &omeg = OwnOpScalar<ConstantCoefficient>(-pow2(omega) * epsilon);
      auto &restr_muinv = OwnOpScalar<RestrictedCoefficient>(muinv, attr);
      auto &restr_omeg = OwnOpScalar<RestrictedCoefficient>(omeg, attr);

      a.AddDomainIntegrator(new CurlCurlIntegrator(restr_muinv), nullptr);
      a.AddDomainIntegrator(new VectorFEMassIntegrator(restr_omeg), nullptr);

      const int cdim = (dim == 2) ? 1 : dim;
      auto &pml_c1_re =
         OwnOpVector<PmlDiagMatrixCoefficient>(cdim, this,
                                               &Ex25PmlContext::DetJInvJtJRe);
      auto &pml_c1_im =
         OwnOpVector<PmlDiagMatrixCoefficient>(cdim, this,
                                               &Ex25PmlContext::DetJInvJtJIm);
      auto &c1_re = OwnOpVector<ScalarVectorProductCoefficient>(muinv, pml_c1_re);
      auto &c1_im = OwnOpVector<ScalarVectorProductCoefficient>(muinv, pml_c1_im);
      auto &restr_c1_re = OwnOpVector<VectorRestrictedCoefficient>(c1_re, attr_pml);
      auto &restr_c1_im = OwnOpVector<VectorRestrictedCoefficient>(c1_im, attr_pml);

      auto &pml_c2_re =
         OwnOpVector<PmlDiagMatrixCoefficient>(dim, this,
                                               &Ex25PmlContext::DetJJtJInvRe);
      auto &pml_c2_im =
         OwnOpVector<PmlDiagMatrixCoefficient>(dim, this,
                                               &Ex25PmlContext::DetJJtJInvIm);
      auto &c2_re = OwnOpVector<ScalarVectorProductCoefficient>(omeg, pml_c2_re);
      auto &c2_im = OwnOpVector<ScalarVectorProductCoefficient>(omeg, pml_c2_im);
      auto &restr_c2_re = OwnOpVector<VectorRestrictedCoefficient>(c2_re, attr_pml);
      auto &restr_c2_im = OwnOpVector<VectorRestrictedCoefficient>(c2_im, attr_pml);

      a.AddDomainIntegrator(new CurlCurlIntegrator(restr_c1_re),
                            new CurlCurlIntegrator(restr_c1_im));
      a.AddDomainIntegrator(new VectorFEMassIntegrator(restr_c2_re),
                            new VectorFEMassIntegrator(restr_c2_im));
   }

   void AddPreconditionerIntegrators(BilinearForm &prec)
   {
      Array<int> attr;
      Array<int> attr_pml;
      if (mesh->attributes.Size())
      {
         attr.SetSize(mesh->attributes.Max());
         attr_pml.SetSize(mesh->attributes.Max());
         attr = 0;
         attr[0] = 1;
         attr_pml = 0;
         if (mesh->attributes.Max() > 1) { attr_pml[1] = 1; }
      }

      auto &muinv = OwnPrecScalar<ConstantCoefficient>(1.0 / mu);
      auto &absomeg = OwnPrecScalar<ConstantCoefficient>(pow2(omega) * epsilon);
      auto &restr_muinv = OwnPrecScalar<RestrictedCoefficient>(muinv, attr);
      auto &restr_absomeg = OwnPrecScalar<RestrictedCoefficient>(absomeg, attr);
      prec.AddDomainIntegrator(new CurlCurlIntegrator(restr_muinv));
      prec.AddDomainIntegrator(new VectorFEMassIntegrator(restr_absomeg));

      const int cdim = (dim == 2) ? 1 : dim;
      auto &pml_c1_abs =
         OwnPrecVector<PmlDiagMatrixCoefficient>(cdim, this,
                                                 &Ex25PmlContext::DetJInvJtJAbs);
      auto &c1_abs = OwnPrecVector<ScalarVectorProductCoefficient>(muinv, pml_c1_abs);
      auto &restr_c1_abs =
         OwnPrecVector<VectorRestrictedCoefficient>(c1_abs, attr_pml);

      auto &pml_c2_abs =
         OwnPrecVector<PmlDiagMatrixCoefficient>(dim, this,
                                                 &Ex25PmlContext::DetJJtJInvAbs);
      auto &c2_abs =
         OwnPrecVector<ScalarVectorProductCoefficient>(absomeg, pml_c2_abs);
      auto &restr_c2_abs =
         OwnPrecVector<VectorRestrictedCoefficient>(c2_abs, attr_pml);

      prec.AddDomainIntegrator(new CurlCurlIntegrator(restr_c1_abs));
      prec.AddDomainIntegrator(new VectorFEMassIntegrator(restr_c2_abs));
   }

   void ComputeErrors(const ComplexGridFunction &x,
                      int order,
                      double *rel_re,
                      double *rel_im,
                      double *total) const
   {
      VectorFunctionCoefficient ex_re(dim, [this](const Vector &xx, Vector &yy)
      {
         ExactReal(xx, yy);
      });
      VectorFunctionCoefficient ex_im(dim, [this](const Vector &xx, Vector &yy)
      {
         ExactImag(xx, yy);
      });
      const int order_quad = std::max(2, 2 * order + 1);
      const mfem::IntegrationRule *irs[Geometry::NumGeom];
      for (int i = 0; i < Geometry::NumGeom; ++i)
      {
         irs[i] = &IntRules.Get(i, order_quad);
      }

      const mfem::real_t l2_re =
         x.real().ComputeL2Error(ex_re,
                                 const_cast<const mfem::IntegrationRule **>(irs),
                                 &pml->GetMarkedPMLElements());
      const mfem::real_t l2_im =
         x.imag().ComputeL2Error(ex_im,
                                 const_cast<const mfem::IntegrationRule **>(irs),
                                 &pml->GetMarkedPMLElements());

      ComplexGridFunction zero(const_cast<mfem::FiniteElementSpace *>(x.FESpace()));
      zero = 0.0;
      const mfem::real_t norm_re =
         zero.real().ComputeL2Error(ex_re,
                                    const_cast<const mfem::IntegrationRule **>(irs),
                                    &pml->GetMarkedPMLElements());
      const mfem::real_t norm_im =
         zero.imag().ComputeL2Error(ex_im,
                                    const_cast<const mfem::IntegrationRule **>(irs),
                                    &pml->GetMarkedPMLElements());

      *rel_re = l2_re / norm_re;
      *rel_im = l2_im / norm_im;
      *total = std::sqrt(l2_re * l2_re + l2_im * l2_im);
   }
};

void PmlDiagMatrixCoefficient::Eval(Vector &K, ElementTransformation &T,
                                    const mfem::IntegrationPoint &ip)
{
   mfem::real_t x[3];
   Vector transip(x, 3);
   T.Transform(ip, transip);
   K.SetSize(vdim);
   (context->*function)(transip, K);
}

} // namespace

extern "C" {

   void *CMFEM_Ex25_NewPml(CMFEM_Mesh *mesh, int problem)
   {
      return new Ex25PmlContext(cmfem::As<mfem::Mesh>(mesh), problem);
   }

   void CMFEM_Ex25_DeletePml(void *context)
   {
      delete static_cast<Ex25PmlContext *>(context);
   }

   void CMFEM_Ex25_SetAttributesPml(void *context, CMFEM_Mesh *mesh)
   {
      static_cast<Ex25PmlContext *>(context)->SetAttributes(
         cmfem::As<mfem::Mesh>(mesh));
   }

   void CMFEM_Ex25_BuildEssentialBoundaryPml(void *context,
                                             const CMFEM_Mesh *mesh,
                                             CMFEM_ArrayInt *ess_bdr)
   {
      static_cast<Ex25PmlContext *>(context)->BuildEssentialBoundary(
         *cmfem::As<const mfem::Mesh>(mesh),
         cmfem::ArrayIntRef(ess_bdr));
   }

   void CMFEM_Ex25_ProjectBoundaryPml(void *context,
                                      double mu,
                                      double epsilon,
                                      double omega,
                                      CMFEM_ComplexGridFunction *x,
                                      const CMFEM_ArrayInt *ess_bdr)
   {
      auto *ctx = static_cast<Ex25PmlContext *>(context);
      ctx->SetParameters(static_cast<mfem::real_t>(mu),
                         static_cast<mfem::real_t>(epsilon),
                         static_cast<mfem::real_t>(omega));
      ctx->ProjectBoundary(*cmfem::As<mfem::ComplexGridFunction>(x),
                           const_cast<mfem::Array<int> &>(cmfem::ArrayIntRef(ess_bdr)));
   }

   void CMFEM_Ex25_AssembleRhsPml(void *context,
                                  double mu,
                                  double epsilon,
                                  double omega,
                                  CMFEM_ComplexLinearForm *b)
   {
      auto *ctx = static_cast<Ex25PmlContext *>(context);
      ctx->SetParameters(static_cast<mfem::real_t>(mu),
                         static_cast<mfem::real_t>(epsilon),
                         static_cast<mfem::real_t>(omega));
      ctx->AssembleRhs(*cmfem::As<mfem::ComplexLinearForm>(b));
   }

   void CMFEM_Ex25_AddOperatorIntegratorsPml(void *context,
                                             double mu,
                                             double epsilon,
                                             double omega,
                                             CMFEM_SesquilinearForm *a)
   {
      auto *ctx = static_cast<Ex25PmlContext *>(context);
      ctx->SetParameters(static_cast<mfem::real_t>(mu),
                         static_cast<mfem::real_t>(epsilon),
                         static_cast<mfem::real_t>(omega));
      ctx->AddOperatorIntegrators(*cmfem::As<mfem::SesquilinearForm>(a));
   }

   void CMFEM_Ex25_AddPreconditionerIntegratorsPml(void *context,
                                                   double mu,
                                                   double epsilon,
                                                   double omega,
                                                   CMFEM_BilinearForm *prec)
   {
      auto *ctx = static_cast<Ex25PmlContext *>(context);
      ctx->SetParameters(static_cast<mfem::real_t>(mu),
                         static_cast<mfem::real_t>(epsilon),
                         static_cast<mfem::real_t>(omega));
      ctx->AddPreconditionerIntegrators(*cmfem::As<mfem::BilinearForm>(prec));
   }

   void CMFEM_Ex25_ComputeErrorsPml(void *context,
                                    double mu,
                                    double epsilon,
                                    double omega,
                                    const CMFEM_ComplexGridFunction *x,
                                    int order,
                                    double *rel_re,
                                    double *rel_im,
                                    double *total)
   {
      auto *ctx = static_cast<Ex25PmlContext *>(context);
      ctx->SetParameters(static_cast<mfem::real_t>(mu),
                         static_cast<mfem::real_t>(epsilon),
                         static_cast<mfem::real_t>(omega));
      ctx->ComputeErrors(*cmfem::As<const mfem::ComplexGridFunction>(x),
                         order,
                         rel_re,
                         rel_im,
                         total);
   }

} // extern "C"
