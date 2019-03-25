// Copyright (c) 2019 Anyar, Inc.
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//      http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace linear
{
   enum class BoundaryMode : int
   {
      Error,
      Constant,
      Linear
   };

   struct BoundaryPolicy
   {
      BoundaryMode lower = BoundaryMode::Error;
      BoundaryMode upper = BoundaryMode::Error;

      void set_all(const BoundaryMode& mode) noexcept
      {
         lower = mode;
         upper = mode;
      }
   };

   template <class value_t, size_t N>
   struct Bounds
   {
      std::array<size_t, N> mid{};
      std::array<value_t, N> weight{};
      bool error = false;
   };

   template <class bounds_t, class value_t, class axis_t>
   void search_axis(bounds_t& bounds, const size_t i, const BoundaryPolicy& policy, const axis_t& axis, const value_t& value) noexcept
   {
      if (std::empty(axis))
         return;

      const size_t n = axis.size();

      if (value >= axis.back())
      {
         switch (policy.upper)
         {
         case BoundaryMode::Error:
            bounds.error = true;
            break;
         case BoundaryMode::Constant:
            bounds.mid[i] = n - 2;
            bounds.weight[i] = 1.0;
            break;
         case BoundaryMode::Linear:
            bounds.mid[i] = n - 2;
            if (axis[n - 1] == 0)
            {
               bounds.weight[i] = 1.0 + std::abs(value / axis[n - 2]);
            }
            else
            {
               bounds.weight[i] = 1.0 + std::abs((value - axis[n - 1]) / axis[n - 1]);
            }
            break;
         default:
            break;
         }
      }
      else if (value < axis.front())
      {
         switch (policy.lower)
         {
         case BoundaryMode::Error:
            bounds.error = true;
            break;
         case BoundaryMode::Constant:
            bounds.mid[i] = 0;
            bounds.weight[i] = 0.0;
            break;
         case BoundaryMode::Linear:
            bounds.mid[i] = 0;
            if (axis[0] == 0)
            {
               bounds.weight[i] = 1.0 + std::abs((value - axis[0]) / axis[1]);
            }
            else
            {
               bounds.weight[i] = 1.0 + std::abs((value - axis[0]) / axis[0]);
            }
            break;
         default:
            break;
         }
      }
      else if (value > axis.front())
      {
         // binary search to find tick
         size_t l = 0, h = n - 2;
         auto& mid = bounds.mid[i];
         auto& weight = bounds.weight[i];
         while (l <= h)
         {
            mid = l + (h - l) / 2;
            if (value < axis[mid]) {
               h = mid - 1;
            }
            else if (value >= axis[mid + 1]) {
               l = mid + 1;
            }
            else {
               weight = (value - axis[mid]) / (axis[mid + 1] - axis[mid]);
               break;
            }
         }
      }
   }

   template <size_t Dimensions>
   struct natord
   {
      size_t operator()(const size_t* nd, const size_t* indices) const noexcept
      {
         size_t i = Dimensions - 1;
         size_t product = 1;
         size_t index = indices[i];
         while (i > 0) {
            product *= nd[i--];
            index += indices[i] * product;
         }
         return index;
      }
   };

   template <size_t Dimensions>
   struct rnatord
   {
      size_t operator()(const size_t* nd, const size_t* indices) const noexcept
      {
         size_t i = 0;
         size_t product = 1;
         size_t index = indices[i];
         while (i < Dimensions - 1) {
            product *= nd[i++];
            index += indices[i] * product;
         }
         return index;
      }
   };

   template <class value_t, size_t N>
   struct LUT
   {
      using value_type = value_t;

      BoundaryPolicy boundary_policy;

      std::vector<value_t> data;

      void add_axis(const std::vector<value_t>& axis)
      {
         axes.emplace_back(axis);
         axes_sizes.emplace_back(axis.size());
      }

      const std::vector<value_t>& axis(const size_t i) const
      {
         return axes[i];
      }

      value_t operator()(const value_t& input) const
      {
         Bounds<value_t, N> bounds;
         search_axis(bounds, 0, boundary_policy, axes[0], input);

         if (bounds.error)
            throw std::runtime_error("out of bounds access");

         const size_t l = bounds.mid[0];
         if (l + 1 < axes[0].size())
            return data[l] + (data[l + 1] - data[l]) * bounds.weight[0];
         else
            return data[l];
      }

      value_t operator()(const std::vector<value_t>& input) const
      {
         if (input.size() != N)
            throw std::runtime_error("mismatching input length for LUT lookup");

         Bounds<value_t, N> bounds;
         for (size_t i = 0; i < N; ++i)
         {
            search_axis(bounds, i, boundary_policy, axes[i], input[i]);
         }

         if (bounds.error)
            throw std::runtime_error((std::string("out of bounds access for LUT<") + typeid(value_t).name() + "," + std::to_string(N) + ">").c_str());

         value_t output{};
         value_t factor;
         std::array<size_t, N> buffer;
         static constexpr size_t Power = static_cast<size_t>(1) << N; // 2^N
         for (size_t s = 0; s < Power; ++s)
         {
            factor = static_cast<value_t>(1);
            for (size_t i = 0; i < N; ++i)
            {
               if (s & (static_cast<size_t>(1) << i)) // check if bit i (starting from 0) is set in s
               {
                  buffer[i] = bounds.mid[i];
                  factor *= 1 - bounds.weight[i];
               }
               else
               {
                  buffer[i] = bounds.mid[i] + 1;
                  factor *= bounds.weight[i];
               }
            }
            if (factor > std::numeric_limits<value_t>::epsilon())
            {
               const size_t k = natord<N>()(axes_sizes.data(), buffer.data());
               output += factor * data[k];
            }
         }

         return output;
      }

      std::vector<std::vector<value_t>> axes;
      std::vector<size_t> axes_sizes;
   };
}