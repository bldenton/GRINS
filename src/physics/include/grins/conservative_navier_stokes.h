//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef GRINS_CONSERVATIVE_NAVIER_STOKES_H
#define GRINS_CONSERVATIVE_NAVIER_STOKES_H

// GRINS
#include "grins/physics.h"
#include "grins/assembly_context.h"
#include "grins/grins_enums.h"
#include "grins/single_variable.h"
#include "grins/multi_component_vector_variable.h"

// libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/point.h"

namespace GRINS
{
  // Physics calls for Conservative Navier Stokes
  /*
      This Physics class implements the conservative form of the Navier-Stokes equations
  */
  
  template<class Viscosity, class SpecificHeat, class ThermalConductivity>
  class ConservativeNavierStokes : public Physics
  {
  public:
    ConservativeNavierStokes(const PhysicsName& physics_name,
                             const GetPot& input );
    
    ~ConservativeNavierStokes(){};
    
    // -----------------------------------------------------------------
    // Overwrite GRINS PHYSICS FUNCTIONS - Polymorhpic
    // -----------------------------------------------------------------
    
    // Sets variables to be time-evolving
    // in this case, density, x-, y-, z-momentum, and conservative energy
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system);
    
    // Context Initialization 
    virtual void init_context( AssemblyContext& context );
    
    
    // Register all parameters in this physics and in its property classes
    virtual void register_parameter( const std::string & param_name,
                                     libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
                                     const;
    
    // Auxiliary Variable Initialization
    virtual void auxiliary_init( MultiphysicsSystem& system );
    
    // Time Dependent Part(s)
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context);
    
    // F(U) calculations                                      
    virtual void assemble_mass_time_derivative( bool compute_jacobian,
                                                AssemblyContext & context);
    
    virtual void assemble_momentum_energy_time_derivative( bool compute_jacobian,
                                                           AssemblyContext & context);
    
    
    // Mass Matrix Part(s)
    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );
    
    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value);
    
    // Register Postprocessing variables
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real>& postprocessing );
    
    
    
    // -----------------------------------------------------------------
    // Define CompressibleNavierStokes Functions
    // -----------------------------------------------------------------
    
    
  protected:
    // Solution Variables
    DensityFEVariable& _density_var;
    ConservativeMomentumVariable& _momentum_vars;
    ConservativeEnergyFEVariable& _conserv_energy_var;
  
    // Viscosity Object
    Viscosity _mu;
    
    // Specific Heat Object
    SpecificHeat _cp;
    
    // Thermal Conductivity Object
    ThermalConductivity _k;
    
    // Ratio of Specific Heat
    libMesh::Real _gamma, _R;
    
    // Gravity Vector
    libMesh::Point _g;
    
    // Initial Density (This is a hack. probably need to remove this in the end)
    libMesh::Real _rho_initial;
    
  private:
    ConservativeNavierStokes();
    
    // Read Options from GetPoot input file
    void read_input_options( const GetPot& input );
    
  };
  
  // ------------------------------------------------------------------
  // Define inline functions for CompressibleNavierStokes Class
  // ------------------------------------------------------------------




}  // namespace GRINS


#endif // GRINS_CONSERVATIVE_NAVIER_STOKES_H