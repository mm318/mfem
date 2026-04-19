#include "common.hpp"

#include <cmath>

namespace
{

void QuadTrans(mfem::real_t u, mfem::real_t v, mfem::real_t radius,
               mfem::real_t &x, mfem::real_t &y)
{
   const mfem::real_t d = 4.0 * radius * (M_SQRT2 - 2.0 * radius) *
                          (1.0 - 2.0 * v);

   const mfem::real_t v0 =
      (1.0 + M_SQRT2) * (M_SQRT2 * radius - 2.0 * v) *
      ((4.0 - 3.0 * M_SQRT2) * radius +
       (8.0 * (M_SQRT2 - 1.0) * radius - 2.0) * v) / d;

   const mfem::real_t r =
      2.0 * ((M_SQRT2 - 1.0) * radius * radius * (1.0 - 4.0 * v) +
             2.0 * (1.0 + M_SQRT2 *
                    (1.0 + 2.0 * (2.0 * radius - M_SQRT2 - 1.0) * radius)) *
             v * v) / d;

   const mfem::real_t t = std::asin(v / r) * u / v;
   x = r * std::sin(t);
   y = r * std::cos(t) - v0;
}

void TransformPeriodicHolePoint(const mfem::Vector &u, mfem::real_t radius,
                                mfem::Vector &x)
{
   const mfem::real_t tol = 1e-4;
   x.SetSize(u.Size());

   if (u[1] > 0.5 - tol || u[1] < -0.5 + tol)
   {
      x = u;
      return;
   }
   if (u[0] > 1.0 - tol || u[0] < -1.0 + tol || std::abs(u[0]) < tol)
   {
      x = u;
      return;
   }

   if (u[0] > 0.0)
   {
      if (u[1] > std::abs(u[0] - 0.5))
      {
         QuadTrans(u[0] - 0.5, u[1], radius, x[0], x[1]);
         x[0] += 0.5;
         return;
      }
      if (u[1] < -std::abs(u[0] - 0.5))
      {
         QuadTrans(u[0] - 0.5, -u[1], radius, x[0], x[1]);
         x[0] += 0.5;
         x[1] *= -1.0;
         return;
      }
      if (u[0] - 0.5 > std::abs(u[1]))
      {
         QuadTrans(u[1], u[0] - 0.5, radius, x[1], x[0]);
         x[0] += 0.5;
         return;
      }
      if (u[0] - 0.5 < -std::abs(u[1]))
      {
         QuadTrans(u[1], 0.5 - u[0], radius, x[1], x[0]);
         x[0] *= -1.0;
         x[0] += 0.5;
         return;
      }
   }
   else
   {
      if (u[1] > std::abs(u[0] + 0.5))
      {
         QuadTrans(u[0] + 0.5, u[1], radius, x[0], x[1]);
         x[0] -= 0.5;
         return;
      }
      if (u[1] < -std::abs(u[0] + 0.5))
      {
         QuadTrans(u[0] + 0.5, -u[1], radius, x[0], x[1]);
         x[0] -= 0.5;
         x[1] *= -1.0;
         return;
      }
      if (u[0] + 0.5 > std::abs(u[1]))
      {
         QuadTrans(u[1], u[0] + 0.5, radius, x[1], x[0]);
         x[0] -= 0.5;
         return;
      }
      if (u[0] + 0.5 < -std::abs(u[1]))
      {
         QuadTrans(u[1], -0.5 - u[0], radius, x[1], x[0]);
         x[0] *= -1.0;
         x[0] -= 0.5;
         return;
      }
   }

   x = u;
}

} // namespace

extern "C" {

   CMFEM_Mesh *CMFEM_NewPeriodicTwoHoleMesh2d(int ref_levels, double radius)
   {
      mfem::Mesh *mesh = new mfem::Mesh(2, 29, 16, 24, 2);
      int vi[4];

      for (int i = 0; i < 2; i++)
      {
         const int o = 13 * i;
         vi[0] = o + 0; vi[1] = o + 3; vi[2] = o + 4; vi[3] = o + 1;
         mesh->AddQuad(vi);

         vi[0] = o + 1; vi[1] = o + 4; vi[2] = o + 5; vi[3] = o + 2;
         mesh->AddQuad(vi);

         vi[0] = o + 5; vi[1] = o + 8; vi[2] = o + 9; vi[3] = o + 2;
         mesh->AddQuad(vi);

         vi[0] = o + 8; vi[1] = o + 12; vi[2] = o + 15; vi[3] = o + 9;
         mesh->AddQuad(vi);

         vi[0] = o + 11; vi[1] = o + 14; vi[2] = o + 15; vi[3] = o + 12;
         mesh->AddQuad(vi);

         vi[0] = o + 10; vi[1] = o + 13; vi[2] = o + 14; vi[3] = o + 11;
         mesh->AddQuad(vi);

         vi[0] = o + 6; vi[1] = o + 13; vi[2] = o + 10; vi[3] = o + 7;
         mesh->AddQuad(vi);

         vi[0] = o + 0; vi[1] = o + 6; vi[2] = o + 7; vi[3] = o + 3;
         mesh->AddQuad(vi);
      }

      vi[0] = 0; vi[1] = 6; mesh->AddBdrSegment(vi, 1);
      vi[0] = 6; vi[1] = 13; mesh->AddBdrSegment(vi, 1);
      vi[0] = 13; vi[1] = 19; mesh->AddBdrSegment(vi, 1);
      vi[0] = 19; vi[1] = 26; mesh->AddBdrSegment(vi, 1);

      vi[0] = 28; vi[1] = 22; mesh->AddBdrSegment(vi, 2);
      vi[0] = 22; vi[1] = 15; mesh->AddBdrSegment(vi, 2);
      vi[0] = 15; vi[1] = 9; mesh->AddBdrSegment(vi, 2);
      vi[0] = 9; vi[1] = 2; mesh->AddBdrSegment(vi, 2);

      for (int i = 0; i < 2; i++)
      {
         const int o = 13 * i;
         vi[0] = o + 7; vi[1] = o + 3; mesh->AddBdrSegment(vi, 3 + i);
         vi[0] = o + 10; vi[1] = o + 7; mesh->AddBdrSegment(vi, 3 + i);
         vi[0] = o + 11; vi[1] = o + 10; mesh->AddBdrSegment(vi, 3 + i);
         vi[0] = o + 12; vi[1] = o + 11; mesh->AddBdrSegment(vi, 3 + i);
         vi[0] = o + 8; vi[1] = o + 12; mesh->AddBdrSegment(vi, 3 + i);
         vi[0] = o + 5; vi[1] = o + 8; mesh->AddBdrSegment(vi, 3 + i);
         vi[0] = o + 4; vi[1] = o + 5; mesh->AddBdrSegment(vi, 3 + i);
         vi[0] = o + 3; vi[1] = o + 4; mesh->AddBdrSegment(vi, 3 + i);
      }

      mfem::real_t d[2];
      const mfem::real_t a = static_cast<mfem::real_t>(radius / M_SQRT2);

      d[0] = -1.0; d[1] = -0.5; mesh->AddVertex(d);
      d[0] = -1.0; d[1] = 0.0; mesh->AddVertex(d);
      d[0] = -1.0; d[1] = 0.5; mesh->AddVertex(d);

      d[0] = -0.5 - a; d[1] = -a; mesh->AddVertex(d);
      d[0] = -0.5 - a; d[1] = 0.0; mesh->AddVertex(d);
      d[0] = -0.5 - a; d[1] = a; mesh->AddVertex(d);

      d[0] = -0.5; d[1] = -0.5; mesh->AddVertex(d);
      d[0] = -0.5; d[1] = -a; mesh->AddVertex(d);
      d[0] = -0.5; d[1] = a; mesh->AddVertex(d);
      d[0] = -0.5; d[1] = 0.5; mesh->AddVertex(d);

      d[0] = -0.5 + a; d[1] = -a; mesh->AddVertex(d);
      d[0] = -0.5 + a; d[1] = 0.0; mesh->AddVertex(d);
      d[0] = -0.5 + a; d[1] = a; mesh->AddVertex(d);

      d[0] = 0.0; d[1] = -0.5; mesh->AddVertex(d);
      d[0] = 0.0; d[1] = 0.0; mesh->AddVertex(d);
      d[0] = 0.0; d[1] = 0.5; mesh->AddVertex(d);

      d[0] = 0.5 - a; d[1] = -a; mesh->AddVertex(d);
      d[0] = 0.5 - a; d[1] = 0.0; mesh->AddVertex(d);
      d[0] = 0.5 - a; d[1] = a; mesh->AddVertex(d);

      d[0] = 0.5; d[1] = -0.5; mesh->AddVertex(d);
      d[0] = 0.5; d[1] = -a; mesh->AddVertex(d);
      d[0] = 0.5; d[1] = a; mesh->AddVertex(d);
      d[0] = 0.5; d[1] = 0.5; mesh->AddVertex(d);

      d[0] = 0.5 + a; d[1] = -a; mesh->AddVertex(d);
      d[0] = 0.5 + a; d[1] = 0.0; mesh->AddVertex(d);
      d[0] = 0.5 + a; d[1] = a; mesh->AddVertex(d);

      d[0] = 1.0; d[1] = -0.5; mesh->AddVertex(d);
      d[0] = 1.0; d[1] = 0.0; mesh->AddVertex(d);
      d[0] = 1.0; d[1] = 0.5; mesh->AddVertex(d);

      mesh->FinalizeTopology();
      mesh->SetCurvature(1, true);

      {
         mfem::Array<int> v2v(mesh->GetNV());
         for (int i = 0; i < v2v.Size() - 3; i++)
         {
            v2v[i] = i;
         }
         v2v[v2v.Size() - 3] = 0;
         v2v[v2v.Size() - 2] = 1;
         v2v[v2v.Size() - 1] = 2;

         for (int i = 0; i < mesh->GetNE(); i++)
         {
            mfem::Element *el = mesh->GetElement(i);
            int *v = el->GetVertices();
            const int nv = el->GetNVertices();
            for (int j = 0; j < nv; j++)
            {
               v[j] = v2v[v[j]];
            }
         }
         for (int i = 0; i < mesh->GetNBE(); i++)
         {
            mfem::Element *el = mesh->GetBdrElement(i);
            int *v = el->GetVertices();
            const int nv = el->GetNVertices();
            for (int j = 0; j < nv; j++)
            {
               v[j] = v2v[v[j]];
            }
         }
         mesh->RemoveUnusedVertices();
         mesh->RemoveInternalBoundaries();
      }

      mesh->SetCurvature(3, true);

      for (int level = 0; level < ref_levels; level++)
      {
         mesh->UniformRefinement();
      }

      mesh->Transform([radius_real = static_cast<mfem::real_t>(radius)](
                         const mfem::Vector &u, mfem::Vector &x)
      {
         TransformPeriodicHolePoint(u, radius_real, x);
      });

      return reinterpret_cast<CMFEM_Mesh *>(mesh);
   }

} // extern "C"
