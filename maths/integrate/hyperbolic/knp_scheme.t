/*
 * knp_scheme.t
 *
 *  Created on: 16 Apr 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

template<int rank, template<int> class Model>
void KurganovNoellePetrova<rank, Model>::setField(int d, Field &field)
{
  fields[d] = &field;
}

template<int rank, template<int> class Model>
inline double KurganovNoellePetrova<rank, Model>::van_leer(double u, double up, double um) const
{
  double du = (up-u)*(u-um);

  return (du>0.0)?du/(up-um):0.0;
}


template<int rank, template<int> class Model>
void KurganovNoellePetrova<rank, Model>::reconstruct(size_t direction, const Index &pos, int dir, FluidValues& u) const
{
  Index posp = pos;
  ++posp[direction];
  Index posm = pos;
  --posm[direction];

  for (size_t d=0; d<dim; ++d)
  {
    Field &F = *fields[d];
    u[d] = F[pos] + dir*van_leer(F[pos], F[posp], F[posm]);
  }
}

template<int rank, template<int> class Model>
void KurganovNoellePetrova<rank, Model>::minmax_local_speed(
        size_t direction,
        const FluidValues &uW,
        const FluidValues &uE,
        const InternalVars &pW,
        const InternalVars &pE,
        double &ap,
        double &am) const
{
  double vW, vE;
  double cfW, cfE;

  vW = this->flow_speed(direction, uW, pW);
  vE = this->flow_speed(direction, uE, pE);

  cfW = this->sound_speed(uW, pW);
  cfE = this->sound_speed(uE, pE);

  ap = std::max( (vW+cfW), std::max( (vE+cfE), 0.0 ));
  am = std::min( (vW-cfW), std::min( (vE-cfE), 0.0 ));
}

template<int rank, template<int> class Model>
inline void KurganovNoellePetrova<rank, Model>::flux(size_t direction, const Index &pos, FluidValues& flux) const
{
  FluidValues uW, uE;
  double ap, am;
  FluidValues fE, fW;
  InternalVars pE, pW;

  Index posp = pos;
  ++posp[direction];

  // reconstruct the MHD variables on the cell boundary
  reconstruct(direction, pos,  +1, uE);
  reconstruct(direction, posp, -1, uW);

  // calculate the thermodynamical variables, pressure and temperature
  this->calc_internal_vars(uE, pE);
  this->calc_internal_vars(uW, pW);

  // determine the minimum and maximum local speeds
  this->minmax_local_speed(direction, uW, uE, pW, pE, ap, am);

  // evaluate the flux function
  this->flux_function(direction, uW, pW, fW);
  this->flux_function(direction, uE, pE, fE);

  // assemble everything to calculate the flux
  flux = (ap*fE - am*fW + ap*am*(uW-uE))/(ap-am);
}

template<int rank, template<int> class Model>
inline void KurganovNoellePetrova<rank, Model>::rhs(Index pos, FluidValues& dudt) const
{

  FluidValues sum = 0;
  for (size_t i=0; i<rank; ++i)
  {
    Index posm = pos;
    --posm[i];

    FluidValues fm, fp;
    flux(i, posm, fm);
    flux(i, pos,  fp);

    sum += (fm - fp) / this->getDx()[i];
  }

  dudt = sum;
}
