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

#include "nlohmann/json.hpp"

namespace linear
{
   template <class lut_t>
   inline void to_json(nlohmann::json& j, const lut_t& lut)
   {
      j = nlohmann::json{ {"axes", lut.axes}, {"axes_sizes", lut.axes_sizes}, {"data", lut.data} };
   }

   template <class lut_t>
   inline void from_json(const nlohmann::json& j, lut_t& lut)
   {
      j.at("axes").get_to(lut.axes);
      j.at("axes_sizes").get_to(lut.axes_sizes);
      j.at("data").get_to(lut.data);
   }
}