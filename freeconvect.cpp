#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <math.h>

using namespace std;

double interpolation(double x, double x1, double y1, double x2, double y2) {
  return y1 + (x - x1)*(y2 - y1)/(x2 - x1);
}

pair<double, double> get_rayleigh(double temperatureK) {
  vector<double> temperatures{ 280, 300, 320, 340, 360, 400 };
  vector<double> prandtl{ 0.713, 0.708, 0.703, 0.699, 0.695, 0.689 };
  vector<double> grash_coeff{ 1.815e8, 1.327e8, 0.9942e8, 0.7502e8, 0.5828e8, 0.3653e8 };

  assert(temperatureK > 280 && temperatureK < 400);

  auto lowT_iter = lower_bound(temperatures.begin(), temperatures.end(), temperatureK);
  auto highT_iter = upper_bound(temperatures.begin(), temperatures.end(), temperatureK);
  if (lowT_iter == highT_iter) --lowT_iter;

  int lowT_idx = 0;
  int highT_idx = 0;
  for (int i = 0; i < temperatures.size(); ++i) {
    if (*lowT_iter == temperatures[i]) lowT_idx = i;
    if (*highT_iter == temperatures[i]) highT_idx = i;
  }

  if (*lowT_iter == temperatureK) {
    return pair<double, double>(prandtl[lowT_idx], grash_coeff[lowT_idx]);
  }

  // cout << lowT_idx << " " << *lowT_iter << " " << highT_idx << " " << *highT_iter << endl;

  double prandtl_inter = interpolation(temperatureK, *lowT_iter, prandtl[lowT_idx], *highT_iter, prandtl[highT_idx]);
  double grash_coeff_inter = interpolation(temperatureK, *lowT_iter, grash_coeff[lowT_idx], *highT_iter, grash_coeff[highT_idx]);
  
  cout << "Prandtl: " << prandtl_inter << ", Grashof Coeff: " << grash_coeff_inter << endl;

  return pair<double, double>(prandtl_inter, grash_coeff_inter);
}

int main() {
  double x_in; // in inches
  double T_wall; // in Celsius
  double T_infinity; // in Celsius

  cout << "x in inches: ";
  cin >> x_in;
  cout << "Temperature of wall in Celsius: ";
  cin >> T_wall;
  cout << "Ambient temperature in Celsius: ";
  cin >> T_infinity;

  // convert to m
  double x_m = x_in * 0.0254;

  double T_film = (T_wall+T_infinity) / 2.0;

  pair<double, double> coeffs = get_rayleigh(T_film + 273.0);
  double prandtl = coeffs.first;
  double grashof = coeffs.second * abs(T_wall - T_infinity) * x_m * x_m * x_m;
  double rayleigh_num = prandtl * grashof;
  cout << rayleigh_num << endl;

  double thermal_bndry_layer = 0.0;
  if (rayleigh_num < 1e9) {
    cout << "Laminar\n";
    thermal_bndry_layer = x_m * 3.93 * pow(prandtl, -0.5) * pow((0.952+prandtl), 0.25) * pow(grashof, -0.25);
    cout << "Calculated thermal boundary layer: " << thermal_bndry_layer << " meters\n";
  }
  else {
    cout << "Turbulent\n";
    thermal_bndry_layer = x_m * 0.565 * pow(grashof, -0.1) * pow(prandtl, -0.533333) * pow((1+0.494*pow(prandtl, 0.666666)), 0.1);
    cout << "Calculated thermal boundary layer: " << thermal_bndry_layer << " meters\n";
  }
}