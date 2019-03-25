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

#include "interpolation/Linear.h"
#include "interpolation/Json.h"

int main()
{
   std::vector<double> x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
   std::vector<double> y = { 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 };

   linear::LUT<double, 1> tester;
   tester.boundary_policy.set_all(linear::BoundaryMode::Linear);
   tester.add_axis(x);
   tester.data = y;
   std::cout << "1D: " << tester(2.5) << '\n';
   std::cout << "1D: " << tester(12.0) << '\n';
   std::cout << "1D: " << tester(-5.0) << '\n';

   linear::LUT<double, 2> tester2d;
   tester2d.boundary_policy.set_all(linear::BoundaryMode::Linear);

   //tester2d.add_axis({ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
   //tester2d.add_axis({ 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 });

   tester2d.add_axis({ -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
   tester2d.add_axis({ -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 });

   //tester2d.add_axis({ -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0 });
   //tester2d.add_axis({ -10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0 });

   for (size_t i = 0; i < tester2d.axis(0).size(); ++i)
   {
      for (size_t j = 0; j < tester2d.axis(1).size(); ++j)
         tester2d.data.emplace_back(tester2d.axis(0)[i] * tester2d.axis(1)[j]);
   }
   std::cout << "2D: " << tester2d({ 1.0, 0.0 }) << '\n';
   std::cout << "2D: " << tester2d({ 2.0, 5.5 }) << '\n';
   std::cout << "2D: " << tester2d({ 2.0, 12.0 }) << '\n';
   std::cout << "2D: " << tester2d({ 5.0, 10.0 }) << '\n';
   std::cout << "2D: " << tester2d({ 5.1, 10.0 }) << '\n';
   std::cout << "2D: " << tester2d({ -2.0, 2.0 }) << '\n';
   std::cout << "2D: " << tester2d({ -2.0, -3.0 }) << '\n';

   nlohmann::json j = tester2d;

   linear::LUT<double, 2> lut = j;
   std::cout << lut({ 3.0, 6.0 }) << '\n';

   std::cout << j << '\n';

   return 0;
}