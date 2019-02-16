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

// This Class
#include "grins/conservative_navier_stokes.h"

// GRINS
#include "grins_config.h"
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/constant_viscosity.h"
#include "grins/constant_specific_heat.h"
#include "grins/constant_conductivity.h"
#include "grins/grins_enums.h"
#include "grins/materials_parsing.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/generic_ic_handler.h"
#include "grins/postprocessed_quantities.h"
#include "grins/physics_naming.h"
#include "grins/inc_nav_stokes_macro.h"
#include "grins/parsed_viscosity.h"
#include "grins/parsed_conductivity.h"


// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"
#include "libmesh/dense_vector.h"

// c++ Write Debugging to file
#include <iostream>
#include <fstream>

namespace GRINS
{
  template<class Mu, class SH, class TC>
  ConservativeNavierStokes<Mu, SH, TC>::ConservativeNavierStokes(const std::string& physics_name,
                                                                 const GetPot& input )
    : Physics(physics_name, input),
      _density_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<DensityFEVariable>(VariablesParsing::density_variable_name(input,
                                                                                                                                     PhysicsNaming::conservative_navier_stokes(),
                                                                                                                                     VariablesParsing::PHYSICS))),
      _momentum_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<ConservativeMomentumVariable>(VariablesParsing::conserv_momentum_variable_name(input,
                                                                                                                                     PhysicsNaming::conservative_navier_stokes(), 
                                                                                                                                     VariablesParsing::PHYSICS))),
      _conserv_energy_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<ConservativeEnergyFEVariable>(VariablesParsing::conserv_energy_variable_name(input,
                                                                                                                                     PhysicsNaming::conservative_navier_stokes(), 
                                                                                                                                     VariablesParsing::PHYSICS))),
      _mu(input, MaterialsParsing::material_name(input, PhysicsNaming::conservative_navier_stokes())),
      _cp(input, MaterialsParsing::material_name(input, PhysicsNaming::conservative_navier_stokes())),
      _k(input, MaterialsParsing::material_name(input, PhysicsNaming::conservative_navier_stokes())) //,
      //_gamma(input, MaterialsParsing::material_name(input, PhysicsNaming::conservative_navier_stokes()))
    {
    
    // ---------------------------------------------------------------------------------
    // ??? add additional intialization function calls here ???
    // ---------------------------------------------------------------------------------
    
    // Handler
    this -> _ic_handler = new GenericICHandler( physics_name, input );

    // Problem Dimensionality Error Message
    if ( this -> _momentum_vars.dim() < 2)
      libmesh_error_msg("ERROR: Conservative Navier Stokes is only valid for 2- or 3D Problems. \
      Make sure you have at least two components in your conservative momentum type variable.");
        
    // Read Input Options from input file
    this -> read_input_options(input);
    
    // Check subdomains for variable consistency
    this -> check_var_subdomain_consistency(_density_var);
    this -> check_var_subdomain_consistency(_momentum_vars);
    this -> check_var_subdomain_consistency(_conserv_energy_var);
    }
    
    // ---------------------------------------------------------------------------------
    // Add Coding for any calculation, reading of input values, etc below
    // ---------------------------------------------------------------------------------
    
    // --- Read Input Options
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu,SH,TC>::read_input_options( const GetPot& input )
    {
      //Read Thermodynamic State Information
      MaterialsParsing::read_property( input,
                                       "GasConstant",
                                       PhysicsNaming::conservative_navier_stokes(),
                                       (*this),
                                       _R);
      MaterialsParsing::read_property( input,
                                       "Gamma",
                                       PhysicsNaming::conservative_navier_stokes(),
                                       (*this),
                                       _gamma);
      MaterialsParsing::read_property( input,
                                       "Density",
                                       PhysicsNaming::conservative_navier_stokes(),
                                       (*this),
                                       _rho_initial);     
      
      //Read Gravity Vector
      unsigned int g_dim = input.vector_variable_size("Physics/"+PhysicsNaming::conservative_navier_stokes()+"/g");
      _g(0) = input("Physics/"+PhysicsNaming::conservative_navier_stokes()+"/g", 0.0, 0);
      
      if( g_dim > 1)
        _g(1) = input("Physics/"+PhysicsNaming::conservative_navier_stokes()+"/g", 0.0, 1);
      
      if( g_dim == 3)
        _g(2) = input("Physics/"+PhysicsNaming::conservative_navier_stokes()+"/g", 0.0, 2);
    }
    
    // --- Auxiliary Inputs Initialization
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::auxiliary_init( MultiphysicsSystem& system )
    {
      // Add Auxiliary Input Initialization code here if necessary
    }
    
    // --- Register Postprocessing Variables (Deals with ouput variables. Possibly where we can postprocess momentum, Pressure, Temperature, etc.
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::register_postprocessing_vars( const GetPot& input,
                                                                             PostProcessedQuantities<libMesh::Real>& postprocessing)
    {
      // ---------------------------------------
      // ------        IMPORTANT          ------
      // ------ NEEDS TO BE UPDATED LATER ------
      // ---------------------------------------
    }
    
    // --- Set Time Evolving Variables
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::set_time_evolving_vars( libMesh::FEMSystem* system)
    {
      // Get Dimension of problem: 2- or 3D
      const unsigned int dim = system->get_mesh().mesh_dimension();
      
      // Set Density Variable as Time Evolving
      system->time_evolving(_density_var.rho(),1);
      
      // Set Conservative Energy Variable as Time Evolving
      system->time_evolving(_conserv_energy_var.conserv_energy(),1);
      
      // Set Momentum Variables as Time Evolving [Note :: The numbers indicate discretization order of time derivative]
      system->time_evolving(_momentum_vars.rho_u(),1);
      
      if ( dim > 1)
        system->time_evolving(_momentum_vars.rho_v(),1);      // 2D Domains
        
      if ( dim == 3)
        system->time_evolving(_momentum_vars.rho_w(),1);      // 3D Domains
        
      // Add addtional variables here if necessary
      
    }
    
    // --- Initialize Variables
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::init_context( AssemblyContext& context)
    {
      // Prerequest all data needed to build the linear system or
      // evaluate a quantity of interest
      
      // Density
      context.get_element_fe(_density_var.rho())->get_JxW();    // Jacobian times Weight Function
      context.get_element_fe(_density_var.rho())->get_phi();    // Test Function
      context.get_element_fe(_density_var.rho())->get_dphi();   // Gradient of Test Function
      context.get_element_fe(_density_var.rho())->get_xyz();    // Location in Physical Space
      
      // Conservative Energy
      context.get_element_fe(_conserv_energy_var.conserv_energy())->get_JxW();     // Jacobian time Weight Function
      context.get_element_fe(_conserv_energy_var.conserv_energy())->get_phi();     // Test Function
      context.get_element_fe(_conserv_energy_var.conserv_energy())->get_dphi();    // Gradient of Test Function
      context.get_element_fe(_conserv_energy_var.conserv_energy())->get_xyz();     // Location in Physical Space
      
      // Conservative Momentum :: Code Assumes these are the same for all momentum directions
      context.get_element_fe(_momentum_vars.rho_u())->get_JxW();      // Jacobian times Weighting Function
      context.get_element_fe(_momentum_vars.rho_u())->get_phi();      // Test Function
      context.get_element_fe(_momentum_vars.rho_u())->get_dphi();     // Gradient of Test Function
      context.get_element_fe(_momentum_vars.rho_u())->get_xyz();      // Location in Physical Space
      
    }
    
    // --- Element Time Derivative Calculations F(u)
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::element_time_derivative( bool compute_jacobian,
                                                                        AssemblyContext & context )
    {    
      // ----------------------------------------------------
      // ----- Add Code for calculating Time Derivative -----
      // ----- Could be broken into multiple functions  -----
      // ----- This contains F(u) equations             -----
      // ----------------------------------------------------
      
      this -> assemble_mass_time_derivative( compute_jacobian, context);
      this -> assemble_momentum_energy_time_derivative( compute_jacobian, context);
      //this -> assemble_conserv_energy_time_derivative( compute_jacobian, context);
    }
    
    // --- Mass Residual Calculations M(u)
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::mass_residual( bool compute_jacobian,
                                                              AssemblyContext & context)
    {
      // --------------------------------------------------------
      // ------ Add Code for calculating Mass Residual      -----
      // ------ M(u)                                        -----
      // --------------------------------------------------------
      
      // Open Debug File
      std::ofstream mass_r;
      mass_r.open("mass_residual_NaN.txt");
      
      // -----------------------------------------------------------------------
      // Element Jacobian * Quadrature weights for interior integration.
      // -----------------------------------------------------------------------
      const std::vector<libMesh::Real> &JxW_density = context.get_element_fe(this->_density_var.rho())->get_JxW();
      const std::vector<libMesh::Real> &JxW_momentum = context.get_element_fe(this->_momentum_vars.rho_u())->get_JxW();
      const std::vector<libMesh::Real> &JxW_energy = context.get_element_fe(this->_conserv_energy_var.conserv_energy())->get_JxW();

      // --------------------------------------------------------------------
      // Get Shape Functions at interior quadrature points for each variable
      // --------------------------------------------------------------------
        /* --- Density   --- */
      const std::vector<std::vector<libMesh::Real> >& rho_phi = 
        context.get_element_fe(this->_density_var.rho())->get_phi();
        
        /* --- Momentum   --- */
      const std::vector<std::vector<libMesh::Real> >& momentum_phi = 
        context.get_element_fe(this->_momentum_vars.rho_u())->get_phi();
      
        /* --- Conservative Energy  --- */
      const std::vector<std::vector<libMesh::Real> >& conserv_energy_phi = 
        context.get_element_fe(this->_conserv_energy_var.conserv_energy())->get_phi();
      
      // ---------------------------------------------------------------------  
      // Get Shape Function Gradients (in global coordiantes)
      // @ interior quadrature points for each 
      // ---------------------------------------------------------------------
        /* --- Density   --- */
      const std::vector<std::vector<libMesh::RealGradient> >& rho_gradphi = 
        context.get_element_fe(this->_density_var.rho())->get_dphi();
        
        /* --- Momentum   --- */
      const std::vector<std::vector<libMesh::RealGradient> >& momentum_gradphi = 
        context.get_element_fe(this->_momentum_vars.rho_u())->get_dphi();
        
        /* --- Conservative Energy   --- */
      const std::vector<std::vector<libMesh::RealGradient> >& conserv_energy_gradphi = 
        context.get_element_fe(this->_conserv_energy_var.conserv_energy())->get_dphi();
      
      // ---------------------------------------------------------------
      // Get the number of local degress of freedom in each variable
      // ---------------------------------------------------------------
        /* ---   Density   --- */
      const unsigned int n_rho_dofs = context.get_dof_indices(this->_density_var.rho()).size();
      
        /* ---   Momentum   --- */
      const unsigned int n_rho_u_dofs = context.get_dof_indices(this->_momentum_vars.rho_u()).size();
      
        /* ---   Check to see if number of dofs for X-Momentum & Y-Momentum are the same   --- */
      libmesh_assert (n_rho_u_dofs == context.get_dof_indices(this->_momentum_vars.rho_v()).size());
      
        /* ---   If 3D Domain, get Z-Momentum   --- */
      if (this -> _momentum_vars.dim() == 3)
        libmesh_assert (n_rho_u_dofs == context.get_dof_indices(this->_momentum_vars.rho_w()).size());
        
        /* ---   Conservative Energy   --- */
      const unsigned int n_conserv_energy_dofs = context.get_dof_indices(this->_conserv_energy_var.conserv_energy()).size();      
      
      // ----------------------------------------------------- 
      // Define Residual Frho_u [ R_{rho_u} ]
      // -----------------------------------------------------
      libMesh::DenseSubVector<libMesh::Real> &Frho = context.get_elem_residual(this->_density_var.rho());            // R_(rho)
      libMesh::DenseSubVector<libMesh::Real> &Frho_u = context.get_elem_residual(this->_momentum_vars.rho_u());      // R_{rho_u}
      libMesh::DenseSubVector<libMesh::Real> &Frho_v = context.get_elem_residual(this->_momentum_vars.rho_v());      // R_{row_v}
      libMesh::DenseSubVector<libMesh::Real>* Frho_w = NULL;
      libMesh::DenseSubVector<libMesh::Real> &Fconserv_energy = context.get_elem_residual(this->_conserv_energy_var.conserv_energy());      // R_{conserv_energy}   
      
      // -----------------------------------------------------
      // Define Jacobians
      // -----------------------------------------------------
      libMesh::DenseSubMatrix<libMesh::Real> &Mrho_rho = context.get_elem_jacobian(this->_density_var.rho(), this->_density_var.rho()); // R_{rho},{rho}      
      libMesh::DenseSubMatrix<libMesh::Real> &Mrho_u_rho_u = context.get_elem_jacobian(this->_momentum_vars.rho_u(), this->_momentum_vars.rho_u()); // R_{rho_u},{rho_u}        
      libMesh::DenseSubMatrix<libMesh::Real> &Mrho_v_rho_v = context.get_elem_jacobian(this->_momentum_vars.rho_v(), this->_momentum_vars.rho_v()); // R_{rho_v},{rho_v}       
      libMesh::DenseSubMatrix<libMesh::Real>* Mrho_w_rho_w = NULL;     
      libMesh::DenseSubMatrix<libMesh::Real> &Menergy_energy = context.get_elem_jacobian(this->_conserv_energy_var.conserv_energy(), this->_conserv_energy_var.conserv_energy()); // R_{conserv_energy},{conserv_energy} 
      
      if (this -> _momentum_vars.dim() == 3)
        {
          Frho_w = &context.get_elem_residual(this->_momentum_vars.rho_w());      // R_{row_w}
          Mrho_w_rho_w = &context.get_elem_jacobian(this->_momentum_vars.rho_w(), this->_momentum_vars.rho_w()); // R_{rho_w},{rho_w}
        }
 
      // ---------------------------------------------------------------------------------------------------------
      // Define Vectors and Matrices to be used for intermediate calculations of the Residual and Jacobians
      // ---------------------------------------------------------------------------------------------------------
      libMesh::DenseVector<libMesh::Number> dUdx(5, 0.), dUdy(5, 0.), dUdz(5, 0.), dWdx(5, 0.), dWdy(5, 0.), dWdz(5, 0.);
      libMesh::DenseVector<libMesh::Number> a1_urow(5, 0.), a1_vrow(5, 0.), a1_wrow(5, 0.), a1_energyrow(5, 0.),
                                            a2_urow(5, 0.), a2_vrow(5, 0.), a2_wrow(5, 0.), a2_energyrow(5, 0.),
                                            a3_urow(5, 0.), a3_vrow(5, 0.), a3_wrow(5, 0.), a3_energyrow(5, 0.);
      libMesh::DenseVector<libMesh::Number> dUdt(5, 0.);
   
      // ---------------------------------------------------------------
      // Get the number of element quadrature points
      // ---------------------------------------------------------------
      unsigned int n_qpoints = context.get_element_qrule().n_points();      

      // ------------------------------------------------------------------------
      // Loop over Quadrature Points and assemble Residuals
      // ------------------------------------------------------------------------      
      for (unsigned int qp = 0; qp != n_qpoints; ++qp)
        {
          // For the mass residual, we need to be a little careful.
          // The time integrator is handling the time-discretization
          // for us so we need to supply M(u_fixed) *u' for the residual.
          // u_fixed will be given by the fixed_interior_value function
          // while u' will be given by the interior_rate functions.
          
          libMesh::Real rho_dot, rho_u_dot, rho_v_dot, rho_w_dot = 0.0, energy_dot;
          
          /* ---  Get Interior time rates for density, etc.  --- */
          context.interior_rate(this -> _density_var.rho(), qp, rho_dot);
          context.interior_rate(this -> _momentum_vars.rho_u(), qp, rho_u_dot);
          context.interior_rate(this -> _momentum_vars.rho_v(), qp, rho_v_dot);
          context.interior_rate(this -> _conserv_energy_var.conserv_energy(), qp, energy_dot);
          
          if (this -> _momentum_vars.dim() == 3)
            context.interior_rate(this -> _momentum_vars.rho_w(), qp, rho_w_dot);

          // --------------------------------------------------
          // Compute the Solutions at the Old Newton Iteration
          // --------------------------------------------------
          libMesh::Number density, u_momentum, v_momentum, w_momentum, conserv_energy;
          
            /* --- Density Value at Quadrature Point --- */
          density = context.interior_value(this->_density_var.rho(), qp);
          
            /* --- Momentum Values at Quadrature Point --- */
          u_momentum = context.interior_value(this->_momentum_vars.rho_u(), qp);
          v_momentum = context.interior_value(this->_momentum_vars.rho_v(), qp);
          if (this->_momentum_vars.dim() == 3)
            w_momentum = context.interior_value(this->_momentum_vars.rho_w(), qp);
            
            /* --- Conservative Energy Value at Quadrature Point --- */
          conserv_energy = context.interior_value(this->_conserv_energy_var.conserv_energy(), qp);
          
          // -----------------------------------------------------------
          // Compute the Solution Gradients at the Old Newton Iteration
          // -----------------------------------------------------------
          libMesh::Gradient grad_density, grad_u_momentum, grad_v_momentum, grad_w_momentum, grad_conserv_energy;
          
            /* --- Density Gradietn at Quadrature Point --- */
          grad_density = context.interior_gradient(this->_density_var.rho(), qp);           // grad(density)
          
            /* --- Momentum Gradient at Quadrature Point --- */
          grad_u_momentum = context.interior_gradient(this->_momentum_vars.rho_u(), qp);    // grad(u_momentum)
          grad_v_momentum = context.interior_gradient(this->_momentum_vars.rho_v(), qp);    // grad(v_momentum)
          if (this->_momentum_vars.dim() == 3)
            grad_w_momentum = context.interior_gradient(this->_momentum_vars.rho_w(), qp);  // grad(w_momentum)
            
            /* --- Conservative Energy at Quadrature Point --- */
          grad_conserv_energy = context.interior_gradient(this->_conserv_energy_var.conserv_energy(), qp);  // grad(conserv_energy)
          
          // -------------------------------------------------------------------
          // Compute Viscoity and Thermal Conductivity at this Quadrature Point
          // -------------------------------------------------------------------
          libMesh::Real _mu_qp = this ->_mu(context, qp);
          libMesh::Real _k_qp = this-> _k(context, qp);
          libMesh::Real _cp_qp = this-> _cp();    //_cp(context, qp);
          libMesh::Real _gamma_qp = _gamma;  //this -> _gamma(context, qp);
          libMesh::Real _R_qp = _R;
          
            /* --- Declare flow aligned length scale variable --- */
          libMesh::Real hvel_qp, stab_SUPG_rho, stab_SUPG_momentum, stab_SUPG_energy;
          
            /* --- Time Step Assumed to be 1.  --- */
          libMesh::Number dtime = 1.;
          
          // ----------------------------------------------------------------------------------------
          // Calculate jacobian ai
          // ----------------------------------------------------------------------------------------
          libMesh::Real inv_density = 1./density;
          libMesh::Real sqr_density = (density * density);
          libMesh::Real sqr_u_momentum = u_momentum * u_momentum;
          libMesh::Real sqr_v_momentum = v_momentum * v_momentum;
          libMesh::Real sqr_w_momentum = (this->_momentum_vars.dim() == 3)?(w_momentum*w_momentum):0;
          
          libMesh::Real lambda = -(2./3.)*_mu_qp;
          libMesh::Real mu_R = 2.*_mu_qp + lambda;
          
            /* --- Calculate Temperature, Pressure and local speed of sound @ quadrature point --- */
          libMesh::Real T_qp = (_gamma_qp/(density*_cp_qp)) * (conserv_energy - ((1./(2.*density)) * (sqr_u_momentum + sqr_v_momentum + sqr_w_momentum)));
          libMesh::Real P_qp = (_gamma_qp - 1.) * (conserv_energy - ((1./(2.*density)) * (sqr_u_momentum + sqr_v_momentum + sqr_w_momentum)));
          libMesh::Real a_qp = pow(_gamma_qp * _R_qp * T_qp, 1./2.);
          
            /* --- Calculate Velocity Vector  ---*/  
          libMesh::Number velocity_vec_length;
          libMesh::DenseVector<libMesh::Number> velocity, unit_velocity;
          velocity_vec_length = pow((1./sqr_density)*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum), 1./2.);
          unit_velocity(0) = u_momentum / velocity_vec_length;
          unit_velocity(1) = v_momentum / velocity_vec_length;
          unit_velocity(2) = w_momentum / velocity_vec_length;            
            
            /* --- calculate a1 matrix  --- */
          a1_urow(0) = (_gamma_qp - 3.) * sqr_u_momentum / (2. * sqr_density) + ((_gamma_qp - 1.)/(2. * sqr_density))*(sqr_v_momentum + sqr_w_momentum);
          a1_urow(1) = (3. - _gamma_qp) * u_momentum/density;
          a1_urow(2) = (1. - _gamma_qp) * v_momentum/density;
          a1_urow(3) = (1. - _gamma_qp) * w_momentum/density;
          a1_urow(4) = (_gamma_qp - 1.);
          
          a1_vrow(0) = -u_momentum * v_momentum / sqr_density;
          a1_vrow(1) = v_momentum / density;
          a1_vrow(2) = u_momentum / density;
          a1_vrow(3) = 0.;
          a1_vrow(4) = 0.;
          
          a1_wrow(0) = -u_momentum * w_momentum / sqr_density;
          a1_wrow(1) = w_momentum / density;
          a1_wrow(2) = 0.;
          a1_wrow(3) = u_momentum / density;
          a1_wrow(4) = 0.;
          
          a1_energyrow(0) = (-_gamma_qp*conserv_energy*u_momentum/sqr_density) + ((_gamma_qp - 1.)/(sqr_density*density))*u_momentum*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum);
          a1_energyrow(1) = (_gamma_qp*conserv_energy/density) + ((1.-_gamma_qp)/(2.*sqr_density))*(3.*sqr_u_momentum + sqr_v_momentum + sqr_w_momentum);
          a1_energyrow(2) = (1.-_gamma_qp)*u_momentum*v_momentum/sqr_density;
          a1_energyrow(3) = (1.-_gamma_qp)*u_momentum*w_momentum/sqr_density;
          a1_energyrow(4) = _gamma_qp*u_momentum/density;
          
            /* --- calculate a2 matrix  --- */
          a2_urow(0) = -u_momentum * v_momentum / sqr_density;
          a2_urow(1) = v_momentum / density;
          a2_urow(2) = u_momentum / density;
          a2_urow(3) = 0.;
          a2_urow(4) = 0.;
          
          a2_vrow(0) = (_gamma_qp - 3.) * sqr_v_momentum / (2. * sqr_density) + ((_gamma_qp - 1.)/(2. * sqr_density))*(sqr_u_momentum * sqr_w_momentum);
          a2_vrow(1) = (1. - _gamma_qp) * u_momentum / density;
          a2_vrow(2) = (3. - _gamma_qp) * v_momentum / density;
          a2_vrow(3) = (1. - _gamma_qp) * w_momentum / density;
          a2_vrow(4) = (_gamma_qp - 1.);
          
          a2_wrow(0) = -v_momentum * w_momentum / sqr_density;
          a2_wrow(1) = 0.;
          a2_wrow(2) = w_momentum / density;
          a2_wrow(3) = v_momentum / density;
          a2_wrow(4) = 0.;
          
          a2_energyrow(0) = (-_gamma_qp*conserv_energy*v_momentum/sqr_density) + ((_gamma_qp - 1.)/(sqr_density*density))*v_momentum*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum);
          a2_energyrow(1) = a1_energyrow(2);
          a2_energyrow(2) = (_gamma_qp*conserv_energy/density) + ((1.-_gamma_qp)/(2.*sqr_density))*(sqr_u_momentum + 3.*sqr_v_momentum + sqr_w_momentum);
          a2_energyrow(3) = (1.-_gamma_qp)*v_momentum*w_momentum/sqr_density;
          a2_energyrow(4) = _gamma_qp*v_momentum/density;
          
            /* ---  calculate a3 matrix  --- */
          a3_urow(0) = -u_momentum * w_momentum / sqr_density;
          a3_urow(1) = w_momentum / density;
          a3_urow(2) = 0.;
          a3_urow(3) = u_momentum / density;
          a3_urow(4) = 0.;
          
          a3_vrow(0) = -v_momentum * w_momentum / sqr_density;
          a3_vrow(1) = 0.;
          a3_vrow(2) = w_momentum / density;
          a3_vrow(3) = v_momentum / density;
          a3_vrow(4) = 0.;
          
          a3_wrow(0) = (_gamma_qp - 3.) * sqr_w_momentum / (2. * sqr_density) + ((_gamma_qp - 1.)/(2. * sqr_density))*(sqr_u_momentum * sqr_v_momentum);
          a3_wrow(1) = (1. - _gamma_qp) * u_momentum / density;
          a3_wrow(2) = (1. - _gamma_qp) * v_momentum / density;
          a3_wrow(3) = (3. - _gamma_qp) * w_momentum / density;
          a3_wrow(4) = (_gamma_qp - 1.); 
          
          a3_energyrow(0) = (-_gamma_qp*conserv_energy*w_momentum/sqr_density) + ((_gamma_qp - 1.)/(sqr_density*density))*w_momentum*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum);
          a3_energyrow(1) = a1_energyrow(3);
          a3_energyrow(2) = a2_energyrow(3);
          a3_energyrow(3) = (_gamma_qp*conserv_energy/density) + ((1.-_gamma_qp)/(2.*sqr_density))*(sqr_u_momentum + sqr_v_momentum + 3.*sqr_w_momentum);
          a3_energyrow(4) = _gamma_qp*w_momentum/density;            
 
         // Set up dUdt vector
         dUdt(0) = rho_dot;
         dUdt(1) = rho_u_dot;
         dUdt(2) = rho_v_dot;
         dUdt(3) = rho_w_dot;
         dUdt(4) = energy_dot;
          
          /* ---  perform integration  --- */
          for (unsigned int ii = 0; ii != n_rho_dofs; ++ii)
            {
               // Calculate flow aligned element length scale
              hvel_qp = 2./(abs(unit_velocity.dot(rho_gradphi[ii][qp]() )));
              
              // Calculate density SUPG Stabilization Factor
              stab_SUPG_rho = pow((pow(2./dtime, 2.) + pow(((2.*(velocity_vec_length + a_qp)) / hvel_qp), 2.)), -1./2.);      // NOTE: assumes dtime = 1. [Steady-State]           
            
              Frho(ii) -= JxW_density[qp] * 
                          // Conservative Navier-Stokes
                          (rho_phi[ii][qp] * rho_dot
                          //SUPG Stabilization
                          + stab_SUPG_rho *
                          (rho_gradphi[ii][qp](0) * rho_u_dot +
                           rho_gradphi[ii][qp](1) * rho_v_dot +
                           rho_gradphi[ii][qp](2) * rho_w_dot)
                          );
              
              if(std::isnan(Frho(ii)))
                { mass_r << "--- mass_residual() ---" << "\n";
                  mass_r << "Frho(ii) = " << Frho(ii) << "\n";
                  mass_r << "ii = " << ii;
                  mass_r << "JxW_density[qp] = " << JxW_density[qp] << "\n";
                  mass_r << "rho_phi[ii][qp] = " << rho_phi[ii][qp] << "\n";
                  mass_r << "rho_dot = " << rho_dot << "\n";
                  mass_r << " -----------------------" << "\n";
                }
              
              // Need to add Mrho_rho calculations here
            }  // end density quadrature loop
          
          for (unsigned int ii = 0; ii != n_rho_u_dofs; ++ii)
            {
              // Calculate Flow Aligned Length Scale
              hvel_qp = 2./(abs(unit_velocity.dot(momentum_gradphi[ii][qp])));
              
              // Calculate SUPG Momentum Stabilization Factor
              stab_SUPG_momentum = pow(pow(2./dtime, 2.) + pow((2.*(velocity_vec_length + a_qp))/hvel_qp, 2.) + pow((4.*_mu_qp)/(density*pow(hvel_qp, 2.)), 2.), -1./2.);      // NOTE: assumes dtime = 1. [Steady-State]
              
              Frho_u(ii) -= JxW_momentum[qp] * 
                            // Conservative Navier-Stokes
                            (momentum_phi[ii][qp] * rho_u_dot
                            // SUPG Stabilization
                            + stab_SUPG_momentum *
                            (momentum_gradphi[ii][qp](0) * a1_urow.dot(dUdt) +
                             momentum_gradphi[ii][qp](1) * a2_urow.dot(dUdt) +
                             momentum_gradphi[ii][qp](2) * a3_urow.dot(dUdt))
                            );
                            
              Frho_v(ii) -= JxW_momentum[qp] * 
                            // Conservative Navier-Stokes
                            (momentum_phi[ii][qp] * rho_v_dot
                             // SUPG Stabilization
                             + stab_SUPG_momentum *
                            (momentum_gradphi[ii][qp](0) * a1_vrow.dot(dUdt) +
                             momentum_gradphi[ii][qp](1) * a2_vrow.dot(dUdt) +
                             momentum_gradphi[ii][qp](2) * a3_vrow.dot(dUdt))                           
                            );
              
              if ( this -> _momentum_vars.dim() == 3)
              (*Frho_w)(ii) -= JxW_momentum[qp] * 
                            // Conservative Navier-Stokes
                            (momentum_phi[ii][qp] * rho_w_dot
                             // SUPG Stabilization
                             + stab_SUPG_momentum *
                            (momentum_gradphi[ii][qp](0) * a1_wrow.dot(dUdt) +
                             momentum_gradphi[ii][qp](1) * a2_wrow.dot(dUdt) +
                             momentum_gradphi[ii][qp](2) * a3_wrow.dot(dUdt))                             
                            );
              
              if(std::isnan(Frho_u(ii)))
                { mass_r << "--- mass_residual() ---" << "\n";
                  mass_r << "Frho_u(ii) = " << Frho_u(ii) << "\n";
                  mass_r << "ii = " << ii;
                  mass_r << "JxW_momentum[qp] = " << JxW_momentum[qp] << "\n";
                  mass_r << "momentum_phi[ii][qp] = " << momentum_phi[ii][qp] << "\n";
                  mass_r << "rho_u_dot = " << rho_u_dot << "\n";
                  mass_r << " -----------------------" << "\n";
                }
              
              if(std::isnan(Frho_v(ii)))
                { mass_r << "--- mass_residual() ---" << "\n";
                  mass_r << "Frho_v(ii) = " << Frho_v(ii) << "\n";
                  mass_r << "ii = " << ii;
                  mass_r << "JxW_momentum[qp] = " << JxW_momentum[qp] << "\n";
                  mass_r << "momentum_phi[ii][qp] = " << momentum_phi[ii][qp] << "\n";
                  mass_r << "rho_v_dot = " << rho_v_dot << "\n";
                  mass_r << " -----------------------" << "\n";
                }
                
              if ( this -> _momentum_vars.dim() == 3)                
                {
                  if(std::isnan((*Frho_w)(ii)))
                    { mass_r << "--- mass_residual() ---" << "\n";
                      mass_r << "Frho_w(ii) = " << (*Frho_w)(ii) << "\n";
                      mass_r << "ii = " << ii;
                      mass_r << "JxW_momentum[qp] = " << JxW_momentum[qp] << "\n";
                      mass_r << "momentum_phi[ii][qp] = " << momentum_phi[ii][qp] << "\n";
                      mass_r << "rho_w_dot = " << rho_w_dot << "\n";
                      mass_r << " -----------------------" << "\n";
                    }
                 } 
              
              // Need to add Mrho_u_rho_u, etc calculations here
            }  // end momentum quadrature loop
            
          for (unsigned int ii = 0; ii != n_conserv_energy_dofs; ++ii)
            {
              // Calculate Flow Aligned Length Scale
              hvel_qp = 2./(abs(unit_velocity.dot(conserv_energy_gradphi[ii][qp])));
              
              // Calculate SUPG Momentum Stabilization Factor
              stab_SUPG_energy = pow(pow(2./dtime, 2.) + pow((2.*(velocity_vec_length + a_qp))/hvel_qp, 2.) + pow((4.*_k_qp)/(density*_cp_qp*pow(hvel_qp, 2.)), 2.), -1./2.);      // NOTE: assumes dtime = 1. [Steady-State]
              
              Fconserv_energy(ii) -= JxW_energy[qp] * 
                        // Conservative Navier-Stokes
                        (conserv_energy_phi[ii][qp] * energy_dot
                        // SUPG Stabilization
                        + stab_SUPG_energy *
                        (conserv_energy_gradphi[ii][qp](0) * a1_energyrow.dot(dUdt) +
                         conserv_energy_gradphi[ii][qp](1) * a2_energyrow.dot(dUdt) +
                         conserv_energy_gradphi[ii][qp](2) * a3_energyrow.dot(dUdt))                        
                        );
              
              if(std::isnan(Fconserv_energy(ii)))
                { mass_r << "--- mass_residual() ---" << "\n";
                  mass_r << "Fconserv_energy(ii) = " << Fconserv_energy(ii) << "\n";
                  mass_r << "ii = " << ii;
                  mass_r << "JxW_energy[qp] = " << JxW_energy[qp] << "\n";
                  mass_r << "conserv_energy_phi[ii][qp] = " << conserv_energy_phi[ii][qp] << "\n";
                  mass_r << "energy_dot = " << energy_dot << "\n";
                  mass_r << " -----------------------" << "\n";
                }
              
              // Need to add Menergy_energy calculations here
            }  // end conservative energy quadrature loop
        }  // end of qp loop
        
        mass_r.close();    
    }
    
    // ----------------------------------------------------------
    // --- Assemble Mass Time Derivative Portion of F(u)
    // ----------------------------------------------------------
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::assemble_mass_time_derivative( bool compute_jacobian,
                                                                              AssemblyContext & context)
    {
      // Element Jacobian * Quadrature weights for interior integration.
      const std::vector<libMesh::Real> &JxW_density = context.get_element_fe(this->_density_var.rho())->get_JxW();
      
      // Get Shape Functions at interior quadrature points for each variable
        /* --- Density   --- */
      const std::vector<std::vector<libMesh::Real> >& rho_phi = 
        context.get_element_fe(this->_density_var.rho())->get_phi();
        
        /* --- Momentum   --- */
      const std::vector<std::vector<libMesh::Real> >& momentum_phi = 
        context.get_element_fe(this->_momentum_vars.rho_u())->get_phi();
       
      // Get Shape Function Gradients (in global coordiantes)
      // @ interior quadrature points for each variable
        /* --- Density   --- */
      const std::vector<std::vector<libMesh::RealGradient> >& rho_gradphi = 
        context.get_element_fe(this->_density_var.rho())->get_dphi();
        
        /* --- Momentum   --- */
      const std::vector<std::vector<libMesh::RealGradient> >& momentum_gradphi = 
        context.get_element_fe(this->_momentum_vars.rho_u())->get_dphi();
      
      // Define Residual Frho [ R_{rho} ]
      libMesh::DenseSubVector<libMesh::Number> &Frho = context.get_elem_residual(this->_density_var.rho());
      
      // Define Jacobians 
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_rho = context.get_elem_jacobian(this->_density_var.rho(), this->_density_var.rho()); // R_{rho},{rho}
      
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_rho_u = context.get_elem_jacobian(this->_density_var.rho(), this->_momentum_vars.rho_u()); // R_{rho},{rho_u}
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_rho_v = context.get_elem_jacobian(this->_density_var.rho(), this->_momentum_vars.rho_v()); // R_{rho},{rho_v}
      libMesh::DenseSubMatrix<libMesh::Number>* Krho_rho_w = NULL;
      
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_conserv_energy = context.get_elem_jacobian(this->_density_var.rho(), this->_conserv_energy_var.conserv_energy()); // R_{rho},{conserv_energy}
      
      if(this -> _momentum_vars.dim() == 3)
        {
          Krho_rho_w = &context.get_elem_jacobian(this->_density_var.rho(), this->_momentum_vars.rho_w()); // R_{rho},{rho_w}
        }  // end of IF momentum dimesion check
      
      // ---------------------------------------------------------------
      // Get the number of element quadrature points
      // ---------------------------------------------------------------
      unsigned int n_qpoints = context.get_element_qrule().n_points();
      
      // ---------------------------------------------------------------
      // Get the number of local degress of freedom in each variable
      // ---------------------------------------------------------------
        /* ---   Density   --- */
      const unsigned int n_rho_dofs = context.get_dof_indices(this->_density_var.rho()).size();
      
        /* ---   Momentum   --- */
      const unsigned int n_rho_u_dofs = context.get_dof_indices(this->_momentum_vars.rho_u()).size();
      
        /* ---   Check to see if number of dofs for X-Momentum & Y-Momentum are the same   --- */
      libmesh_assert (n_rho_u_dofs == context.get_dof_indices(this->_momentum_vars.rho_v()).size());
      
        /* ---   If 3D Domain, get Z-Momentum   --- */
      if (this -> _momentum_vars.dim() == 3)
        libmesh_assert (n_rho_u_dofs == context.get_dof_indices(this->_momentum_vars.rho_w()).size());
        
        /* ---   Conservative Energy   --- */
      const unsigned int n_conserv_energy_dofs = context.get_dof_indices(this->_conserv_energy_var.conserv_energy()).size();
      
      // ------------------------------------------------------------------------
      // Loop over Quadrature Points and assemble Frho [ R_{rho} ]
      // ------------------------------------------------------------------------
      for (unsigned int qp=0; qp != n_qpoints; qp++)
        {
          // -------------------------------------------------------------------
          // Compute Viscoity and Thermal Conductivity at this Quadrature Point
          // -------------------------------------------------------------------
          libMesh::Real _cp_qp = this-> _cp();    //_cp(context, qp);
          libMesh::Real _gamma_qp = _gamma;  //this -> _gamma(context, qp);
          libMesh::Real _R_qp = _R; 
        
          // --------------------------------------------------------------------
          // Compute the Solution & its gradients at the Old Newton Iteration
          // --------------------------------------------------------------------
            /* ---  Solution  --- */
          libMesh::Number density, u_momentum, v_momentum, w_momentum, conserv_energy;
          density = context.interior_value(this->_density_var.rho(), qp);
          u_momentum = context.interior_value(this->_momentum_vars.rho_u(), qp);
          v_momentum = context.interior_value(this->_momentum_vars.rho_v(), qp);
          if (this->_momentum_vars.dim() == 3)
            w_momentum = context.interior_value(this->_momentum_vars.rho_w(), qp);
          conserv_energy = context.interior_value(this->_conserv_energy_var.conserv_energy(), qp);          
          
            /* --- Calculate Velocity Vector  ---*/  
          libMesh::Number velocity_vec_length, sqr_density, sqr_u_momentum, sqr_v_momentum, sqr_w_momentum;
          libMesh::DenseVector<libMesh::Number> velocity, unit_velocity;
          sqr_density = density * density;
          sqr_u_momentum = u_momentum * u_momentum;
          sqr_v_momentum = v_momentum * v_momentum;
          sqr_w_momentum = w_momentum * w_momentum;
          velocity_vec_length = pow((1./sqr_density)*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum), 1./2.);
          unit_velocity(0) = u_momentum / velocity_vec_length;
          unit_velocity(1) = v_momentum / velocity_vec_length;
          unit_velocity(2) = w_momentum / velocity_vec_length;
          
            /* --- Calculate Temperature, Pressure and local speed of sound @ quadrature point --- */
          libMesh::Real T_qp = (_gamma_qp/(density*_cp_qp)) * (conserv_energy - ((1./(2.*density)) * (sqr_u_momentum + sqr_v_momentum + sqr_w_momentum)));
          libMesh::Real P_qp = (_gamma_qp - 1.) * (conserv_energy - ((1./(2.*density)) * (sqr_u_momentum + sqr_v_momentum + sqr_w_momentum)));
          libMesh::Real a_qp = pow(_gamma_qp*_R_qp*T_qp, 1./2.);          
          
            /* --- Gradients  --- */
          libMesh::Gradient grad_u_momentum, grad_v_momentum, grad_w_momentum;
          grad_u_momentum = context.interior_gradient(this->_momentum_vars.rho_u(), qp);    // grad(u_momentum)
          grad_v_momentum = context.interior_gradient(this->_momentum_vars.rho_v(), qp);    // grad(v_momentum)
          if (this->_momentum_vars.dim() == 3)
            grad_w_momentum = context.interior_gradient(this->_momentum_vars.rho_w(), qp);  // grad(w_momentum)
            
          const libMesh::Number  grad_u_momentum_x = grad_u_momentum(0);
          const libMesh::Number  grad_v_momentum_y = grad_v_momentum(1);
          const libMesh::Number  grad_w_momentum_z = 0.;                  // If 2D Domain         
          if (this->_momentum_vars.dim() == 3)
            const libMesh::Number  grad_z_momentum_z = grad_w_momentum(2);  // If 3D Domain
          
          // Flow Aligned Length Scale
          libMesh::Number hvel_qp, stab_SUPG_rho;
          
          // Time Step Assumed to be 1.
          libMesh::Number dtime = 1.;
          
          // ------------------------------------------------------------
          // Loop over the density degrees of freedom
          // This computes the contributions of the continuity equation
          // ------------------------------------------------------------
          for (unsigned ii=0; ii != n_rho_dofs; ii++)
            {
              // Calculate flow aligned element length scale
              hvel_qp = 2./(abs(unit_velocity.dot(momentum_gradphi[ii][qp])));
              
              // Calculate density SUPG Stabilization Factor
              stab_SUPG_rho = pow(pow(2./dtime, 2.) + pow((2.*(velocity_vec_length + a_qp))/hvel_qp, 2.), -1./2.);      // NOTE: assumes dtime = 1. [Steady-State]
            
              Frho(ii) -= JxW_density[qp] * 
                          // Conservative Navier-Stokes
                          (rho_phi[ii][qp]*(grad_u_momentum_x + grad_v_momentum_y + grad_w_momentum_z) +
                          // SUPG Stabilization
                          stab_SUPG_rho * (rho_gradphi[ii][qp](0)*grad_u_momentum_x + rho_gradphi[ii][qp](1)*grad_v_momentum_y + rho_gradphi[ii][qp](2)*grad_w_momentum_z)
                          );
              
              if(std::isnan(Frho(ii)))
                { std::cout << "--- assemble_mass_time_derivative() ---" << "\n"
                            << "Frho(ii) = " << Frho(ii) << "\n"
                            << "ii = " << ii << "\n"
                            << "JxW_density[qp] = " << JxW_density[qp] << "\n"
                            << "rho_phi[ii][qp] = " << rho_phi[ii][qp] << "\n"
                            << "grad_u_momentum_x = " << grad_u_momentum_x << "\n"
                            << "grad_v_momentum_y = " << grad_v_momentum_y << "\n"
                            << "grad_w_momentum_z = " << grad_w_momentum_z << "\n"
                            << "stabSUPG = " << stab_SUPG_rho << "\n"
                            << " -----------------------" << "\n";
                }
              
              /* if (compute_jacobian)
              {
                // --------------------------------------------
                // Update R_{rho},{rho} Jacobian
                // --------------------------------------------
                for (unsigned int jj=0; jj != n_rho_dofs; jj++)
                  {
                    Krho_rho(ii, jj) -= JxW[qp] * 0.;
                  }  // End of inner jj Loop for R_{rho},{rho} Jacobian Update
                
                // -------------------------------------------------------------
                // Update R_{rho},{rho_u}, R_{rho],{rho_v} and R_{rho},{rho_w}
                // -------------------------------------------------------------  
                for (unsigned int jj=0; jj != n_rho_u_dofs; jj++)
                  {
                    Krho_rho_u(ii, jj) -= JxW[qp] * rho_phi[ii][qp] * 1. * momentum_gradphi[jj][qp](0);
                    Krho_rho_v(ii, jj) -= JxW[qp] * rho_phi[ii][qp] * 1. * momentum_gradphi[jj][qp](1);
                    
                    if (this->_momentum_vars.dim() == 3)
                      {
                        Krho_rho_w(ii, jj) -= JxW[qp] * rho_phi[ii][qp] * 1. * momentum_gradphi[jj][qp](2);
                      }
                  }  // End of inner jj Loop for R_{rho},{rho_u}, R_{rho},{rho_v} and R_{rho},{rho_w}
                
                // ---------------------------------------
                // Update R_{rho},{conserv_energy}
                // ---------------------------------------
                for (unsigned int jj=0; jj != n_conserv_energy_dofs; jj++)
                  {
                    Krho_conserv_energy(ii, jj) -= JxW[qp] * 0.;
                  }  // End of inner jj Loop for R_{rho},{conserv_energy}
              }    // End of compute_jacobian if statment
              */
            }    // End of density degree of freedom loop                      
        }    // End Quadrature Loop  
      
      return;     
    }
    
    // -----------------------------------------------------------------
    // --- Assemble Momentum & Energy Time Derivative Portion of F(u)
    // -----------------------------------------------------------------
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::assemble_momentum_energy_time_derivative( bool compute_jacobian,
                                                                                         AssemblyContext & context)
    {    
      // -----------------------------------------------------------------------
      // Element Jacobian * Quadrature weights for interior integration.
      // -----------------------------------------------------------------------
      const std::vector<libMesh::Real> &JxW_momentum = context.get_element_fe(this->_momentum_vars.rho_u())->get_JxW();
      const std::vector<libMesh::Real> &JxW_energy = context.get_element_fe(this->_conserv_energy_var.conserv_energy())->get_JxW();
      
      // --------------------------------------------------------------------
      // Get Shape Functions at interior quadrature points for each variable
      // --------------------------------------------------------------------
        /* --- Density   --- */
      const std::vector<std::vector<libMesh::Real> >& rho_phi = 
        context.get_element_fe(this->_density_var.rho())->get_phi();
        
        /* --- Momentum   --- */
      const std::vector<std::vector<libMesh::Real> >& momentum_phi = 
        context.get_element_fe(this->_momentum_vars.rho_u())->get_phi();
      
        /* --- Conservative Energy  --- */
      const std::vector<std::vector<libMesh::Real> >& conserv_energy_phi = 
        context.get_element_fe(this->_conserv_energy_var.conserv_energy())->get_phi();
        
      // ------------------------------------------------------------- 
      // Get Shape Function Gradients (in global coordiantes)
      // @ interior quadrature points for each variable
      // -------------------------------------------------------------
        /* --- Density   --- */
      const std::vector<std::vector<libMesh::RealGradient> >& rho_gradphi = 
        context.get_element_fe(this->_density_var.rho())->get_dphi();
        
        /* --- Momentum   --- */
      const std::vector<std::vector<libMesh::RealGradient> >& momentum_gradphi = 
        context.get_element_fe(this->_momentum_vars.rho_u())->get_dphi();
        
        /* --- Conservative Energy  --- */
      const std::vector<std::vector<libMesh::RealGradient> >& conserv_energy_gradphi = 
        context.get_element_fe(this->_conserv_energy_var.conserv_energy())->get_dphi();
      
      // ----------------------------------------------------- 
      // Define Residual Frho_u [ R_{rho_u} ]
      // -----------------------------------------------------
      libMesh::DenseSubVector<libMesh::Number> &Frho_u = context.get_elem_residual(this->_momentum_vars.rho_u());      // R_{rho_u}
      libMesh::DenseSubVector<libMesh::Number> &Frho_v = context.get_elem_residual(this->_momentum_vars.rho_v());      // R_{row_v}
      libMesh::DenseSubVector<libMesh::Number>* Frho_w = NULL;
      libMesh::DenseSubVector<libMesh::Number> &Fconserv_energy = context.get_elem_residual(this->_conserv_energy_var.conserv_energy());      // R_{conserv_energy}
      
      // -----------------------------------------------------
      // Define Jacobians
      // -----------------------------------------------------
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_u_rho = context.get_elem_jacobian(this->_momentum_vars.rho_u(), this->_density_var.rho()); // R_{rho_u},{rho} 
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_u_rho_u = context.get_elem_jacobian(this->_momentum_vars.rho_u(), this->_momentum_vars.rho_u()); // R_{rho_u},{rho_u}     
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_u_rho_v = context.get_elem_jacobian(this->_momentum_vars.rho_u(), this->_momentum_vars.rho_v()); // R_{rho_u},{rho_v}
      libMesh::DenseSubMatrix<libMesh::Number>* Krho_u_rho_w = NULL;
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_u_conserv_energy = context.get_elem_jacobian(this->_momentum_vars.rho_u(), this->_conserv_energy_var.conserv_energy()); // R_{rho_u},{conserv_energy}

      libMesh::DenseSubMatrix<libMesh::Number> &Krho_v_rho = context.get_elem_jacobian(this->_momentum_vars.rho_v(), this->_density_var.rho()); // R_{rho_v},{rho}       
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_v_rho_u = context.get_elem_jacobian(this->_momentum_vars.rho_v(), this->_momentum_vars.rho_u()); // R_{rho_v},{rho_u}     
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_v_rho_v = context.get_elem_jacobian(this->_momentum_vars.rho_v(), this->_momentum_vars.rho_v()); // R_{rho_v},{rho_v}
      libMesh::DenseSubMatrix<libMesh::Number>* Krho_v_rho_w = NULL; 
      libMesh::DenseSubMatrix<libMesh::Number> &Krho_v_conserv_energy = context.get_elem_jacobian(this->_momentum_vars.rho_v(), this->_conserv_energy_var.conserv_energy()); // R_{rho_v},{conserv_energy}

      libMesh::DenseSubMatrix<libMesh::Number>* Krho_w_rho = NULL;
      libMesh::DenseSubMatrix<libMesh::Number>* Krho_w_rho_u = NULL;     
      libMesh::DenseSubMatrix<libMesh::Number>* Krho_w_rho_v = NULL;
      libMesh::DenseSubMatrix<libMesh::Number>* Krho_w_rho_w = NULL;
      libMesh::DenseSubMatrix<libMesh::Number>* Krho_w_conserv_energy = NULL;
      
      libMesh::DenseSubMatrix<libMesh::Number> &Kconserv_energy_rho = context.get_elem_jacobian(this->_conserv_energy_var.conserv_energy(), this->_density_var.rho()); // R_{conserv_energy},{rho} 
      libMesh::DenseSubMatrix<libMesh::Number> &Kconserv_energy_rho_u = context.get_elem_jacobian(this->_conserv_energy_var.conserv_energy(), this->_momentum_vars.rho_u()); // R_{conserv_energy},{rho_u}     
      libMesh::DenseSubMatrix<libMesh::Number> &Kconserv_energy_rho_v = context.get_elem_jacobian(this->_conserv_energy_var.conserv_energy(), this->_momentum_vars.rho_v()); // R_{conserv_energy},{rho_v}
      libMesh::DenseSubMatrix<libMesh::Number>* Kconserv_energy_rho_w = NULL;
      libMesh::DenseSubMatrix<libMesh::Number> &Kconserv_energy_conserv_energy = context.get_elem_jacobian(this->_conserv_energy_var.conserv_energy(), this->_conserv_energy_var.conserv_energy()); // R_{conserv_energy},{conserv_energy}     
      
      if(this -> _momentum_vars.dim() == 3)
        {
          Frho_w  = &context.get_elem_residual(this->_momentum_vars.rho_w());       // R_{rho_w}
        
          Krho_u_rho_w = &context.get_elem_jacobian(this->_momentum_vars.rho_u(), this->_momentum_vars.rho_w()); // R_{rho_u},{rho_w}
          
          Krho_v_rho_w = &context.get_elem_jacobian(this->_momentum_vars.rho_v(), this->_momentum_vars.rho_w()); // R_{rho_v},{rho_w}
          
          Krho_w_rho = &context.get_elem_jacobian(this->_momentum_vars.rho_w(), this->_density_var.rho()); // R_{rho_w},{rho}       
          Krho_w_rho_u = &context.get_elem_jacobian(this->_momentum_vars.rho_w(), this->_momentum_vars.rho_u()); // R_{rho_w},{rho_u}     
          Krho_w_rho_v = &context.get_elem_jacobian(this->_momentum_vars.rho_w(), this->_momentum_vars.rho_v()); // R_{rho_w},{rho_v}
          Krho_w_rho_w = &context.get_elem_jacobian(this->_momentum_vars.rho_w(), this->_momentum_vars.rho_w()); // R_{rho_w},{rho_w}         
          Krho_w_conserv_energy = &context.get_elem_jacobian(this->_momentum_vars.rho_w(), this->_conserv_energy_var.conserv_energy()); // R_{rho_w},{conserv_energy}
          
          Kconserv_energy_rho_w = &context.get_elem_jacobian(this->_conserv_energy_var.conserv_energy(), this->_momentum_vars.rho_w()); // R_{conserv_energy},{rho_w}
        }  // end of IF momentum dimesion == 3 check
      
      // -----------------------------------------------  
      // Get the number of element quadrature points
      // -----------------------------------------------
      unsigned int n_qpoints = context.get_element_qrule().n_points();
      
      // --------------------------------------------------------------
      // Get the number of local degress of freedom in each variable
      // --------------------------------------------------------------
        /* ---   Density   --- */
      const unsigned int n_rho_dofs = context.get_dof_indices(this->_density_var.rho()).size();
      
        /* ---   Momentum   --- */
      const unsigned int n_rho_u_dofs = context.get_dof_indices(this->_momentum_vars.rho_u()).size();
      
        /* ---   Check to see if number of dofs for X-Momentum & Y-Momentum are the same   --- */
      libmesh_assert (n_rho_u_dofs == context.get_dof_indices(this->_momentum_vars.rho_v()).size());
      
        /* ---   If 3D Domain, get Z-Momentum   --- */
      if (this -> _momentum_vars.dim() == 3)
        libmesh_assert (n_rho_u_dofs == context.get_dof_indices(this->_momentum_vars.rho_w()).size());
        
        /* ---   Conservative Energy   --- */
      const unsigned int n_conserv_energy_dofs = context.get_dof_indices(this->_conserv_energy_var.conserv_energy()).size();    
      
      // ---------------------------------------------------------------------------------------------------------
      // Define Vectors and Matrices to be used for intermediate calculations of the Residual and Jacobians
      // ---------------------------------------------------------------------------------------------------------
      libMesh::DenseVector<libMesh::Number> dUdx(5, 0.), dUdy(5, 0.), dUdz(5, 0.), dWdx(5, 0.), dWdy(5, 0.), dWdz(5, 0.);
      libMesh::DenseVector<libMesh::Number> a1_urow(5, 0.), a1_vrow(5, 0.), a1_wrow(5, 0.), a1_energyrow(5, 0.),
                                            a2_urow(5, 0.), a2_vrow(5, 0.), a2_wrow(5, 0.), a2_energyrow(5, 0.),
                                            a3_urow(5, 0.), a3_vrow(5, 0.), a3_wrow(5, 0.), a3_energyrow(5, 0.),
                                            b1_urow(5, 0.), b1_vrow(5, 0.), b1_wrow(5, 0.), b1_energyrow(5, 0.),
                                            b2_urow(5, 0.), b2_vrow(5, 0.), b2_wrow(5, 0.), b2_energyrow(5, 0.),
                                            b3_urow(5, 0.), b3_vrow(5, 0.), b3_wrow(5, 0.), b3_energyrow(5, 0.),
                                            c11_urow(5, 0.), c11_vrow(5, 0.), c11_wrow(5, 0.), c11_energyrow(5, 0.),
                                            c12_urow(5, 0.), c12_vrow(5, 0.), c12_wrow(5, 0.), c12_energyrow(5, 0.),
                                            c13_urow(5, 0.), c13_vrow(5, 0.), c13_wrow(5, 0.), c13_energyrow(5, 0.),
                                            c21_urow(5, 0.), c21_vrow(5, 0.), c21_wrow(5, 0.), c21_energyrow(5, 0.),
                                            c22_urow(5, 0.), c22_vrow(5, 0.), c22_wrow(5, 0.), c22_energyrow(5, 0.),
                                            c23_urow(5, 0.), c23_vrow(5, 0.), c23_wrow(5, 0.), c23_energyrow(5, 0.),
                                            c31_urow(5, 0.), c31_vrow(5, 0.), c31_wrow(5, 0.), c31_energyrow(5, 0.),
                                            c32_urow(5, 0.), c32_vrow(5, 0.), c32_wrow(5, 0.), c32_energyrow(5, 0.),
                                            c33_urow(5, 0.), c33_vrow(5, 0.), c33_wrow(5, 0.), c33_energyrow(5, 0.);
      libMesh::DenseVector<libMesh::Number> d_dx_rho(5, 0.), d_dx_umomentum(5, 0.), d_dx_vmomentum(5, 0.), d_dx_wmomentum(5, 0.), d_dx_conserv_energy(5, 0.),
                                            d_dy_rho(5, 0.), d_dy_umomentum(5, 0.), d_dy_vmomentum(5, 0.), d_dy_wmomentum(5, 0.), d_dy_conserv_energy(5, 0.),
                                            d_dz_rho(5, 0.), d_dz_umomentum(5, 0.), d_dz_vmomentum(5, 0.), d_dz_wmomentum(5, 0.), d_dz_conserv_energy(5, 0.);
      libMesh::DenseVector<libMesh::Number> NS_stab_u(5, 0.), NS_stab_v(5, 0.), NS_stab_w(5, 0.), NS_stab_energy(5, 0.);
      
      // -------------------------------------------------------------------------------------------------------------------------------------------------           
      // Loop over Quadrature Points and assemble
      // Frho [ R_{rho} ], Frho_u [ R_{rho_u} ], Frho_v [R_{rho_v} ], Frho_w [R_{rho_w} ],
      // Krho_u_rho [ R_{rho_u},{rho} ], Krho_u_rho_u [ R_{rho_u},{rho_u} ], Krho_u_rho_v [ R_{rho_u},{rho_v} ], Krho_u_rho_w [ R_{rho_u},{rho_w} ], etc.
      // -------------------------------------------------------------------------------------------------------------------------------------------------
      for (unsigned int qp=0; qp != n_qpoints; qp++)
        {
          // --------------------------------------------------
          // Compute the Solutions at the Old Newton Iteration
          // --------------------------------------------------
          libMesh::Number density, u_momentum, v_momentum, w_momentum, conserv_energy;
          
            /* --- Density Value at Quadrature Point --- */
          density = context.interior_value(this->_density_var.rho(), qp);
          
            /* --- Momentum Values at Quadrature Point --- */
          u_momentum = context.interior_value(this->_momentum_vars.rho_u(), qp);
          v_momentum = context.interior_value(this->_momentum_vars.rho_v(), qp);
          if (this->_momentum_vars.dim() == 3)
            w_momentum = context.interior_value(this->_momentum_vars.rho_w(), qp);
            
            /* --- Conservative Energy Value at Quadrature Point --- */
          conserv_energy = context.interior_value(this->_conserv_energy_var.conserv_energy(), qp);
          
          // -----------------------------------------------------------
          // Compute the Solution Gradients at the Old Newton Iteration
          // -----------------------------------------------------------
          libMesh::Gradient grad_density, grad_u_momentum, grad_v_momentum, grad_w_momentum, grad_conserv_energy;
          
            /* --- Density Gradietn at Quadrature Point --- */
          grad_density = context.interior_gradient(this->_density_var.rho(), qp);           // grad(density)
          
            /* --- Momentum Gradient at Quadrature Point --- */
          grad_u_momentum = context.interior_gradient(this->_momentum_vars.rho_u(), qp);    // grad(u_momentum)
          grad_v_momentum = context.interior_gradient(this->_momentum_vars.rho_v(), qp);    // grad(v_momentum)
          if (this->_momentum_vars.dim() == 3)
            grad_w_momentum = context.interior_gradient(this->_momentum_vars.rho_w(), qp);  // grad(w_momentum)
            
            /* --- Conservative Energy at Quadrature Point --- */
          grad_conserv_energy = context.interior_gradient(this->_conserv_energy_var.conserv_energy(), qp);  // grad(conserv_energy)
          
          // -------------------------------------------------------------------
          // Compute Viscoity and Thermal Conductivity at this Quadrature Point
          // -------------------------------------------------------------------
          libMesh::Real _mu_qp = this ->_mu(context, qp);
          libMesh::Real _k_qp = this-> _k(context, qp);
          libMesh::Real _cp_qp = this-> _cp();    //_cp(context, qp);
          libMesh::Real _gamma_qp = _gamma;  //this -> _gamma(context, qp);
          libMesh::Real _R_qp = _R;
          
            /* --- Declare flow aligned length scale variable --- */
          libMesh::Real hvel_qp, stab_SUPG_momentum, stab_SUPG_energy;
          
            /* --- Time Step Assumed to be 1.  --- */
          libMesh::Number dtime = 1.;
          
          // ----------------------------------------------------------------------------------------
          // Calculate jacobians (ai, bi & cij) Note: Be Conscious of negative sign when using cij
          // ----------------------------------------------------------------------------------------
          libMesh::Real inv_density = 1./density;
          libMesh::Real sqr_density = (density * density);
          libMesh::Real sqr_u_momentum = u_momentum * u_momentum;
          libMesh::Real sqr_v_momentum = v_momentum * v_momentum;
          libMesh::Real sqr_w_momentum = (this->_momentum_vars.dim() == 3)?(w_momentum*w_momentum):0;
          
          libMesh::Real lambda = -(2./3.)*_mu_qp;
          libMesh::Real mu_R = 2.*_mu_qp + lambda;
          
            /* --- Calculate Temperature, Pressure and local speed of sound @ quadrature point --- */
          libMesh::Real T_qp = (_gamma_qp/(density*_cp_qp)) * (conserv_energy - ((1./(2.*density)) * (sqr_u_momentum + sqr_v_momentum + sqr_w_momentum)));
          libMesh::Real P_qp = (_gamma_qp - 1.) * (conserv_energy - ((1./(2.*density)) * (sqr_u_momentum + sqr_v_momentum + sqr_w_momentum)));
          libMesh::Real a_qp = pow(_gamma_qp*_R_qp*T_qp, 1./2.);
          
            /* --- Calculate Velocity Vector  ---*/  
          libMesh::Number velocity_vec_length;
          libMesh::DenseVector<libMesh::Number> velocity, unit_velocity;
          velocity_vec_length = pow((1./sqr_density)*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum), 1./2.);
          unit_velocity(0) = u_momentum / velocity_vec_length;
          unit_velocity(1) = v_momentum / velocity_vec_length;
          unit_velocity(2) = w_momentum / velocity_vec_length;
          
            /* --- calculate viscosity tensor --- */
          libMesh::Real tau_11, tau_12, tau_13, tau_21, tau_22, tau_23, tau_31, tau_32, tau_33;
          
          tau_11 = _mu_qp * ((4./3.)*(inv_density*grad_u_momentum(0) - (u_momentum/sqr_density)*grad_density(0)) - 
                             (2./3.)*(inv_density*grad_v_momentum(1) - (v_momentum/sqr_density)*grad_density(1) + 
                                      inv_density*grad_w_momentum(2) - (w_momentum/sqr_density)*grad_density(2)));
          tau_12 = _mu_qp * (inv_density*grad_u_momentum(1) - (u_momentum/sqr_density)*grad_density(1) + 
                             inv_density*grad_v_momentum(0) - (v_momentum/sqr_density)*grad_density(0));
          tau_13 = _mu_qp * (inv_density*grad_u_momentum(2) - (u_momentum/sqr_density)*grad_density(2) + 
                             inv_density*grad_w_momentum(0) - (w_momentum/sqr_density)*grad_density(0));
          tau_21 = tau_12;
          tau_22 = (2./3.)*_mu_qp * (2.*(inv_density*grad_v_momentum(1) - (v_momentum/sqr_density)*grad_density(1)) - 
                                     inv_density*grad_u_momentum(0) + (u_momentum/sqr_density)*grad_density(0) - 
                                     inv_density*grad_w_momentum(2) + (w_momentum/sqr_density)*grad_density(2));
          tau_23 = _mu_qp * (inv_density*grad_w_momentum(1) - (w_momentum/sqr_density)*grad_density(1) +
                             inv_density*grad_v_momentum(2) - (v_momentum/sqr_density)*grad_density(2));
          tau_31 = tau_13;
          tau_32 = tau_23;
          tau_33 = (2./3.)*_mu_qp * (2.*(inv_density*grad_w_momentum(2) - (w_momentum/sqr_density)*grad_density(2)) - 
                                     inv_density*grad_u_momentum(0) + (u_momentum/sqr_density)*grad_density(0) - 
                                     inv_density*grad_v_momentum(1) + (v_momentum/sqr_density)*grad_density(1));
                   
            /* --- calculate a1 matrix  --- */
          a1_urow(0) = (_gamma_qp - 3.) * sqr_u_momentum / (2. * sqr_density) + ((_gamma_qp - 1.)/(2. * sqr_density))*(sqr_v_momentum + sqr_w_momentum);
          a1_urow(1) = (3. - _gamma_qp) * u_momentum/density;
          a1_urow(2) = (1. - _gamma_qp) * v_momentum/density;
          a1_urow(3) = (1. - _gamma_qp) * w_momentum/density;
          a1_urow(4) = (_gamma_qp - 1.);
          
          a1_vrow(0) = -u_momentum * v_momentum / sqr_density;
          a1_vrow(1) = v_momentum / density;
          a1_vrow(2) = u_momentum / density;
          a1_vrow(3) = 0.;
          a1_vrow(4) = 0.;
          
          a1_wrow(0) = -u_momentum * w_momentum / sqr_density;
          a1_wrow(1) = w_momentum / density;
          a1_wrow(2) = 0.;
          a1_wrow(3) = u_momentum / density;
          a1_wrow(4) = 0.;
          
          a1_energyrow(0) = (-_gamma_qp*conserv_energy*u_momentum/sqr_density) + ((_gamma_qp - 1.)/(sqr_density*density))*u_momentum*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum);
          a1_energyrow(1) = (_gamma_qp*conserv_energy/density) + ((1.-_gamma_qp)/(2.*sqr_density))*(3.*sqr_u_momentum + sqr_v_momentum + sqr_w_momentum);
          a1_energyrow(2) = (1.-_gamma_qp)*u_momentum*v_momentum/sqr_density;
          a1_energyrow(3) = (1.-_gamma_qp)*u_momentum*w_momentum/sqr_density;
          a1_energyrow(4) = _gamma_qp*u_momentum/density;
          
            /* --- calculate a2 matrix  --- */
          a2_urow(0) = -u_momentum * v_momentum / sqr_density;
          a2_urow(1) = v_momentum / density;
          a2_urow(2) = u_momentum / density;
          a2_urow(3) = 0.;
          a2_urow(4) = 0.;
          
          a2_vrow(0) = (_gamma_qp - 3.) * sqr_v_momentum / (2. * sqr_density) + ((_gamma_qp - 1.)/(2. * sqr_density))*(sqr_u_momentum * sqr_w_momentum);
          a2_vrow(1) = (1. - _gamma_qp) * u_momentum / density;
          a2_vrow(2) = (3. - _gamma_qp) * v_momentum / density;
          a2_vrow(3) = (1. - _gamma_qp) * w_momentum / density;
          a2_vrow(4) = (_gamma_qp - 1.);
          
          a2_wrow(0) = -v_momentum * w_momentum / sqr_density;
          a2_wrow(1) = 0.;
          a2_wrow(2) = w_momentum / density;
          a2_wrow(3) = v_momentum / density;
          a2_wrow(4) = 0.;
          
          a2_energyrow(0) = (-_gamma_qp*conserv_energy*v_momentum/sqr_density) + ((_gamma_qp - 1.)/(sqr_density*density))*v_momentum*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum);
          a2_energyrow(1) = a1_energyrow(2);
          a2_energyrow(2) = (_gamma_qp*conserv_energy/density) + ((1.-_gamma_qp)/(2.*sqr_density))*(sqr_u_momentum + 3.*sqr_v_momentum + sqr_w_momentum);
          a2_energyrow(3) = (1.-_gamma_qp)*v_momentum*w_momentum/sqr_density;
          a2_energyrow(4) = _gamma_qp*v_momentum/density;
          
            /* ---  calculate a3 matrix  --- */
          a3_urow(0) = -u_momentum * w_momentum / sqr_density;
          a3_urow(1) = w_momentum / density;
          a3_urow(2) = 0.;
          a3_urow(3) = u_momentum / density;
          a3_urow(4) = 0.;
          
          a3_vrow(0) = -v_momentum * w_momentum / sqr_density;
          a3_vrow(1) = 0.;
          a3_vrow(2) = w_momentum / density;
          a3_vrow(3) = v_momentum / density;
          a3_vrow(4) = 0.;
          
          a3_wrow(0) = (_gamma_qp - 3.) * sqr_w_momentum / (2. * sqr_density) + ((_gamma_qp - 1.)/(2. * sqr_density))*(sqr_u_momentum * sqr_v_momentum);
          a3_wrow(1) = (1. - _gamma_qp) * u_momentum / density;
          a3_wrow(2) = (1. - _gamma_qp) * v_momentum / density;
          a3_wrow(3) = (3. - _gamma_qp) * w_momentum / density;
          a3_wrow(4) = (_gamma_qp - 1.); 
          
          a3_energyrow(0) = (-_gamma_qp*conserv_energy*w_momentum/sqr_density) + ((_gamma_qp - 1.)/(sqr_density*density))*w_momentum*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum);
          a3_energyrow(1) = a1_energyrow(3);
          a3_energyrow(2) = a2_energyrow(3);
          a3_energyrow(3) = (_gamma_qp*conserv_energy/density) + ((1.-_gamma_qp)/(2.*sqr_density))*(sqr_u_momentum + sqr_v_momentum + 3.*sqr_w_momentum);
          a3_energyrow(4) = _gamma_qp*w_momentum/density;
          
            /* --- calculate b1 matrix --- */
          b1_urow(0) = (1./sqr_density) * (mu_R*grad_u_momentum(0) + lambda*(grad_v_momentum(1)+grad_w_momentum(2)) - (mu_R/density)*(2.*u_momentum*grad_density(0) - v_momentum*grad_density(1) - w_momentum*grad_density(2)));
          b1_urow(1) = (mu_R/sqr_density) * grad_density(0);
          b1_urow(2) = (lambda/sqr_density) * grad_density(1);
          b1_urow(3) = (lambda/sqr_density) * grad_density(2);
          b1_urow(4) = 0.;
          
          b1_vrow(0) = (_mu_qp/sqr_density) * (grad_u_momentum(1) + grad_v_momentum(0) - ((2.*u_momentum)/density)*grad_density(1) - ((2.*v_momentum)/density)*grad_density(0));
          b1_vrow(1) = (_mu_qp/sqr_density) * grad_density(1);
          b1_vrow(2) = (_mu_qp/sqr_density) * grad_density(0);
          b1_vrow(3) = 0.;
          b1_vrow(4) = 0.;
          
          b1_wrow(0) = (_mu_qp/sqr_density) * (grad_u_momentum(2) + grad_w_momentum(0) - ((2.*u_momentum)/density)*grad_density(2) - ((2.*w_momentum)/density)*grad_density(0));
          b1_wrow(1) = (_mu_qp/sqr_density) * grad_density(2);
          b1_wrow(2) = 0.;
          b1_wrow(3) = (_mu_qp/sqr_density) * grad_density(0);
          b1_wrow(4) = 0.;
          
          b1_energyrow(0) = inv_density * ((u_momentum/density)*tau_11 + (v_momentum/density)*tau_12 + (w_momentum/density)*tau_13) +
                            (u_momentum/density)*b1_urow(0) + (v_momentum/density)*b1_vrow(0) + (w_momentum/density)*b1_wrow(0)
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * (-grad_conserv_energy(0) + 
                              ((2.*conserv_energy/density) - (3.*sqr_u_momentum/sqr_density) - (3.*sqr_v_momentum/sqr_density) - (3.*sqr_w_momentum/sqr_density))*grad_density(0) +
                              (2.*u_momentum/density)*grad_u_momentum(0) + (2.*v_momentum/density)*grad_v_momentum(0) + (2.*w_momentum/density)*grad_w_momentum(0));
          b1_energyrow(1) = (u_momentum/density)*b1_urow(1) + (v_momentum/density)*b1_vrow(1) + (w_momentum/density)*b1_wrow(1) - tau_11/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*u_momentum/density)*grad_density(0) - grad_u_momentum(0));
          b1_energyrow(2) = (u_momentum/density)*b1_urow(2) + (v_momentum/density)*b1_vrow(2) + (w_momentum/density)*b1_wrow(2) - tau_12/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*v_momentum/density)*grad_density(0) - grad_v_momentum(0));
          b1_energyrow(3) = (u_momentum/density)*b1_urow(3) + (v_momentum/density)*b1_vrow(3) + (w_momentum/density)*b1_wrow(3) - tau_13/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*w_momentum/density)*grad_density(0) - grad_w_momentum(0));
          b1_energyrow(4) = ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * grad_density(0);
          
            /* --- calculate b2 matrix --- */
          b2_urow(0) = b1_vrow(0);
          b2_urow(1) = b1_vrow(1);
          b2_urow(2) = b1_vrow(2);
          b2_urow(3) = 0.;
          b2_urow(4) = 0.;
          
          b2_vrow(0) = (1./sqr_density) * (lambda*(grad_u_momentum(0) + grad_w_momentum(2)) + mu_R*(grad_v_momentum(1) + (u_momentum/density)*grad_density(0) -
                          (2.*v_momentum/density)*grad_density(1) + (w_momentum/density)*grad_density(2)));
          b2_vrow(1) = (lambda/sqr_density)*grad_density(0);
          b2_vrow(2) = (mu_R/sqr_density)*grad_density(1);
          b2_vrow(3) = (lambda/sqr_density)*grad_density(2);
          b2_vrow(4) = 0.;
          
          b2_wrow(0) = (_mu_qp/sqr_density) * (grad_v_momentum(2) + grad_w_momentum(1) - (2.*v_momentum/density)*grad_density(2) - (2.*w_momentum/density)*grad_density(1));
          b2_wrow(1) = 0.;
          b2_wrow(2) = (_mu_qp/sqr_density) * grad_density(2);
          b2_wrow(3) = (_mu_qp/sqr_density) * grad_density(1);
          b2_wrow(4) = 0.;
          
          b2_energyrow(0) = inv_density * ((u_momentum/density)*tau_21 + (v_momentum/density)*tau_22 + (w_momentum/density)*tau_23) +
                            (u_momentum/density)*b2_urow(0) + (v_momentum/density)*b2_vrow(0) + (w_momentum/density)*b2_wrow(0)
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * (-grad_conserv_energy(1) + 
                              ((2.*conserv_energy/density) - (3.*sqr_u_momentum/sqr_density) - (3.*sqr_v_momentum/sqr_density) - (3.*sqr_w_momentum/sqr_density))*grad_density(1) +
                              (2.*u_momentum/density)*grad_u_momentum(1) + (2.*v_momentum/density)*grad_v_momentum(1) + (2.*w_momentum/density)*grad_w_momentum(1));
          b2_energyrow(1) = (u_momentum/density)*b2_urow(1) + (v_momentum/density)*b2_vrow(1) + (w_momentum/density)*b2_wrow(1) - tau_21/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*u_momentum/density)*grad_density(1) - grad_u_momentum(1));
          b2_energyrow(2) = (u_momentum/density)*b2_urow(2) + (v_momentum/density)*b2_vrow(2) + (w_momentum/density)*b2_wrow(2) - tau_22/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*v_momentum/density)*grad_density(1) - grad_v_momentum(1));
          b2_energyrow(3) = (u_momentum/density)*b2_urow(3) + (v_momentum/density)*b2_vrow(3) + (w_momentum/density)*b2_wrow(3) - tau_23/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*w_momentum/density)*grad_density(1) - grad_w_momentum(1));
          b2_energyrow(4) = ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * grad_density(1);
          
            /* --- calculate b3 matrix --- */
          b3_urow(0) = b1_wrow(0);
          b3_urow(1) = b1_wrow(1);
          b3_urow(2) = 0.;
          b3_urow(3) = b1_wrow(3);
          b3_urow(4) = 0.;
          
          b3_vrow(0) = b2_wrow(0);
          b3_vrow(1) = 0.;
          b3_vrow(2) = b2_wrow(2);
          b3_vrow(3) = b2_wrow(3);
          b3_vrow(4) = 0.;
          
          b3_wrow(0) = (1./sqr_density) * (lambda*(grad_u_momentum(0) + grad_v_momentum(1)) + mu_R*(grad_w_momentum(2) + (u_momentum/density)*grad_density(0) -
                          (2.*w_momentum/density)*grad_density(2) + (v_momentum/density)*grad_density(1)));
          b3_wrow(1) = (lambda/sqr_density) * grad_density(0);
          b3_wrow(2) = (lambda/sqr_density) * grad_density(1);
          b3_wrow(3) = (mu_R/sqr_density) * grad_density(2);
          b3_wrow(4) = 0.;
          
          b3_energyrow(0) = inv_density * ((u_momentum/density)*tau_31 + (v_momentum/density)*tau_32 + (w_momentum/density)*tau_33) +
                            (u_momentum/density)*b3_urow(0) + (v_momentum/density)*b3_vrow(0) + (w_momentum/density)*b3_wrow(0)
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * (-grad_conserv_energy(2) + 
                              ((2.*conserv_energy/density) - (3.*sqr_u_momentum/sqr_density) - (3.*sqr_v_momentum/sqr_density) - (3.*sqr_w_momentum/sqr_density))*grad_density(2) +
                              (2.*u_momentum/density)*grad_u_momentum(2) + (2.*v_momentum/density)*grad_v_momentum(2) + (2.*w_momentum/density)*grad_w_momentum(2));
          b3_energyrow(1) = (u_momentum/density)*b3_urow(1) + (v_momentum/density)*b3_vrow(1) + (w_momentum/density)*b3_wrow(1) - tau_31/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*u_momentum/density)*grad_density(2) - grad_u_momentum(2));
          b3_energyrow(2) = (u_momentum/density)*b3_urow(2) + (v_momentum/density)*b3_vrow(2) + (w_momentum/density)*b3_wrow(2) - tau_32/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*v_momentum/density)*grad_density(2) - grad_v_momentum(2));
          b3_energyrow(3) = (u_momentum/density)*b3_urow(3) + (v_momentum/density)*b3_vrow(3) + (w_momentum/density)*b3_wrow(3) - tau_33/density
                            - ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * ((2.*w_momentum/density)*grad_density(2) - grad_w_momentum(2));
          b3_energyrow(4) = ((_k_qp*_cp_qp)/(sqr_density*_gamma_qp)) * grad_density(2);    
          
            /* --- calculate c11 matrix  --- */
          c11_urow(0) = -mu_R*u_momentum/sqr_density;
          c11_urow(1) = mu_R/density;
          c11_urow(2) = 0.;
          c11_urow(3) = 0.;
          c11_urow(4) = 0.;
          
          c11_vrow(0) = -_mu_qp*v_momentum/sqr_density;
          c11_vrow(1) = 0.;
          c11_vrow(2) = _mu_qp/density;
          c11_vrow(3) = 0.;
          c11_vrow(4) = 0.;
          
          c11_wrow(0) = -_mu_qp*w_momentum/sqr_density;
          c11_wrow(1) = 0.;
          c11_wrow(2) = 0.;
          c11_wrow(3) = _mu_qp*density;
          c11_wrow(4) = 0.;
          
          c11_energyrow(0) = (u_momentum/density)*c11_urow(0) + (v_momentum/density)*c11_vrow(0) + (w_momentum/density)*c11_wrow(0) +
                             (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*(-conserv_energy + (1./density)*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum));
          c11_energyrow(1) = (u_momentum/density)*c11_urow(1) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*u_momentum;
          c11_energyrow(2) = (v_momentum/density)*c11_vrow(2) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*v_momentum;
          c11_energyrow(3) = (w_momentum/density)*c11_wrow(3) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*w_momentum;
          c11_energyrow(4) = (_k_qp*_gamma_qp)/(density*_cp_qp);      
            
            /* --- calculate c12 matrix  --- */
          c12_urow(0) = -lambda*v_momentum/sqr_density;
          c12_urow(1) = 0.;
          c12_urow(2) = lambda/density;
          c12_urow(3) = 0.;
          c12_urow(4) = 0.;
          
          c12_vrow(0) = -_mu_qp*u_momentum/sqr_density;
          c12_vrow(1) = _mu_qp/density;
          c12_vrow(2) = 0.;
          c12_vrow(3) = 0.;
          c12_vrow(4) = 0.;
          
          c12_wrow(0) = 0.;
          c12_wrow(1) = 0.;
          c12_wrow(2) = 0.;
          c12_wrow(3) = 0.;
          c12_wrow(4) = 0.;
          
          c12_energyrow(0) = (u_momentum/density)*c12_urow(0) + (v_momentum/density)*c12_vrow(0);
          c12_energyrow(1) = (v_momentum/density)*c12_vrow(1);
          c12_energyrow(2) = (u_momentum/density)*c12_urow(2);
          c12_energyrow(3) = 0.;
          c12_energyrow(4) = 0.;
                    
            /* --- calculate c13 matrix  --- */
          c13_urow(0) = -lambda*w_momentum/sqr_density;
          c13_urow(1) = 0.;
          c13_urow(2) = 0.;
          c13_urow(3) = lambda/density;
          c13_urow(4) = 0.;
          
          c13_vrow(0) = 0.;
          c13_vrow(1) = 0.;
          c13_vrow(2) = 0.;
          c13_vrow(3) = 0.;
          c13_vrow(4) = 0.;
          
          c13_wrow(0) = -_mu_qp*u_momentum/sqr_density;
          c13_wrow(1) = _mu_qp/density;
          c13_wrow(2) = 0.;
          c13_wrow(3) = 0.;
          c13_wrow(4) = 0.;
          
          c13_energyrow(0) = (u_momentum/density)*c13_urow(0) + (w_momentum/density)*c13_wrow(0);
          c13_energyrow(1) = (w_momentum/density)*c13_wrow(1);
          c13_energyrow(2) = 0.;
          c13_energyrow(3) = (u_momentum/density)*c13_urow(3);
          c13_energyrow(4) = 0.;          
            
            /* --- calculate c21 matrix  --- */ 
          c21_urow(0) = c11_vrow(0);
          c21_urow(1) = 0.;
          c21_urow(2) = c11_vrow(2);
          c21_urow(3) = 0.;
          c21_urow(4) = 0.;
          
          c21_vrow(0) = -lambda*u_momentum/sqr_density;
          c21_vrow(1) = lambda/density;
          c21_vrow(2) = 0.;
          c21_vrow(3) = 0.;
          c21_vrow(4) = 0.;
          
          c21_wrow(0) = 0.;
          c21_wrow(1) = 0.;
          c21_wrow(2) = 0.;
          c21_wrow(3) = 0.;
          c21_wrow(4) = 0.;
          
          c21_energyrow(0) = (u_momentum/density)*c21_urow(0) + (v_momentum/density)*c21_vrow(0);
          c21_energyrow(1) = (v_momentum/density)*c21_vrow(1);
          c21_energyrow(2) = (u_momentum/density)*c21_urow(2);
          c21_energyrow(3) = 0.;
          c21_energyrow(4) = 0.;          
          
            /* --- calculate c22 matrix  --- */
          c22_urow(0) = c12_vrow(0);
          c22_urow(1) = c12_vrow(1);
          c22_urow(2) = 0.;
          c22_urow(3) = 0.;
          c22_urow(4) = 0.;
          
          c22_vrow(0) = -mu_R*v_momentum/sqr_density;
          c22_vrow(1) = 0.;
          c22_vrow(2) = mu_R/density;
          c22_vrow(3) = 0.;
          c22_vrow(4) = 0.;
          
          c22_wrow(0) = -_mu_qp*w_momentum/sqr_density;
          c22_wrow(1) = 0.;
          c22_wrow(2) = 0.;
          c22_wrow(3) = _mu_qp/density;
          c22_wrow(4) = 0.;
          
          c22_energyrow(0) = (u_momentum/density)*c22_urow(0) + (v_momentum/density)*c22_vrow(0) + (w_momentum/density)*c22_wrow(0) +
                             (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*(-conserv_energy + (1./density)*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum));
          c22_energyrow(1) = (u_momentum/density)*c22_urow(1) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*u_momentum;
          c22_energyrow(2) = (v_momentum/density)*c22_vrow(2) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*v_momentum;
          c22_energyrow(3) = (w_momentum/density)*c22_wrow(3) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*w_momentum;
          c22_energyrow(4) = (_k_qp*_gamma_qp)/(density*_cp_qp);          
            
            /* --- calculate c23 matrix  --- */
          c23_urow(0) = 0.;
          c23_urow(1) = 0.;
          c23_urow(2) = 0.;
          c23_urow(3) = 0.;
          c23_urow(4) = 0.;
          
          c23_vrow(0) = -lambda*w_momentum/sqr_density;
          c23_vrow(1) = 0.;
          c23_vrow(2) = 0.;
          c23_vrow(3) = lambda/density;
          c23_vrow(4) = 0.;
          
          c23_wrow(0) = -_mu_qp*v_momentum/sqr_density;
          c23_wrow(1) = 0.;
          c23_wrow(2) = _mu_qp/density;
          c23_wrow(3) = 0.;
          c23_wrow(4) = 0.;
          
          c23_energyrow(0) = (v_momentum/density)*c23_vrow(0) + (w_momentum/density)*c23_wrow(0);
          c23_energyrow(1) = 0.;
          c23_energyrow(2) = (w_momentum/density)*c23_wrow(2);
          c23_energyrow(3) = (v_momentum/density)*c23_vrow(3);
          c23_energyrow(4) = 0.;
            
            /* --- calculate c31 matrix  --- */
          c31_urow(0) = c11_wrow(0);
          c31_urow(1) = 0.;
          c31_urow(2) = 0.;
          c31_urow(3) = c11_wrow(3);
          c31_urow(4) = 0.;
          
          c31_vrow(0) = 0.;
          c31_vrow(1) = 0.;
          c31_vrow(2) = 0.;
          c31_vrow(3) = 0.;
          c31_vrow(4) = 0.;
          
          c31_wrow(0) = -lambda*u_momentum/sqr_density;
          c31_wrow(1) = lambda/density;
          c31_wrow(2) = 0.;
          c31_wrow(3) = 0.;
          c31_wrow(4) = 0.;
          
          c31_energyrow(0) = (u_momentum/density)*c31_urow(0) + (w_momentum/density)*c31_wrow(0);
          c31_energyrow(1) = (w_momentum/density)*c31_wrow(1);
          c31_energyrow(2) = 0.;
          c31_energyrow(3) = (u_momentum/density)*c31_urow(3);
          c31_energyrow(4) = 0.;
            
            /* --- calculate c32 matrix  --- */
          c32_urow(0) = 0.;
          c32_urow(1) = 0.;
          c32_urow(2) = 0.;
          c32_urow(3) = 0.;
          c32_urow(4) = 0.;
          
          c32_vrow(0) = c22_wrow(0);
          c32_vrow(1) = 0.;
          c32_vrow(2) = 0.;
          c32_vrow(3) = c22_wrow(3);
          c32_vrow(4) = 0.;
          
          c32_wrow(0) = -lambda*v_momentum/sqr_density;
          c32_wrow(1) = 0.;
          c32_wrow(2) = lambda/density;
          c32_wrow(3) = 0.;
          c32_wrow(4) = 0.;
          
          c32_energyrow(0) = (v_momentum/density)*c32_vrow(0) + (w_momentum/density)*c32_wrow(0);
          c32_energyrow(1) = 0.;
          c32_energyrow(2) = (w_momentum/density)*c32_wrow(2);
          c32_energyrow(3) = (v_momentum/density)*c32_vrow(3);
          c32_energyrow(4) = 0.;          
            
            /* --- calculate c33 matrix  --- */
          c33_urow(0) = c13_wrow(0);
          c33_urow(1) = c13_wrow(1);
          c33_urow(2) = 0.;
          c33_urow(3) = 0.;
          c33_urow(4) = 0.;
          
          c33_vrow(0) = c23_wrow(0);
          c33_vrow(1) = 0.;
          c33_vrow(2) = c23_wrow(2);
          c33_vrow(3) = 0.;
          c33_vrow(4) = 0.;
          
          c33_wrow(0) = -mu_R*w_momentum/sqr_density;
          c33_wrow(1) = 0.;
          c33_wrow(2) = 0.;
          c33_wrow(3) = mu_R/density;
          c33_wrow(4) = 0.; 
          
          c33_energyrow(0) = (u_momentum/density)*c33_urow(0) + (v_momentum/density)*c33_vrow(0) + (w_momentum/density)*c33_wrow(0) +
                             (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*(-conserv_energy + (1./density)*(sqr_u_momentum + sqr_v_momentum + sqr_w_momentum));
          c33_energyrow(1) = (u_momentum/density)*c33_urow(1) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*u_momentum;
          c33_energyrow(2) = (v_momentum/density)*c33_vrow(2) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*v_momentum;
          c33_energyrow(3) = (w_momentum/density)*c33_wrow(3) - (_k_qp*_gamma_qp)/(sqr_density*_cp_qp)*w_momentum;
          c33_energyrow(4) = (_k_qp*_gamma_qp)/(density*_cp_qp);                     
          
          // ------------------------------------------------------------------
          // Organize previous iteration vectors to calculate Residuals
          // U = (rho  rho_u  rho_v  rho_w  conserv_energy)
          // ------------------------------------------------------------------
          dUdx(0) = grad_density(0);
          dUdx(1) = grad_u_momentum(0);
          dUdx(2) = grad_v_momentum(0);
          dUdx(3) = grad_w_momentum(0);
          dUdx(4) = grad_conserv_energy(0);
          
          dUdy(0) = grad_density(1);
          dUdy(1) = grad_u_momentum(1);
          dUdy(2) = grad_v_momentum(1);
          dUdy(3) = grad_w_momentum(1);
          dUdy(4) = grad_conserv_energy(1);
          
          dUdz(0) = grad_density(2);
          dUdz(1) = grad_u_momentum(2);
          dUdz(2) = grad_v_momentum(2);
          dUdz(3) = grad_w_momentum(2);
          dUdz(4) = grad_conserv_energy(2);    
          
          // ------------------------------------------------------------------
          // First, an ii loop over the velocity degrees of freedom. 
          // ------------------------------------------------------------------
          for (unsigned int ii=0; ii != n_rho_u_dofs; ii++)
            {
              // Calculate Flow Aligned Length Scale
              hvel_qp = 2./(abs(unit_velocity.dot(momentum_gradphi[ii][qp])));
              
              // Calculate SUPG Momentum Stabilization Factor
              stab_SUPG_momentum = (pow(2./dtime, 2.) + pow((2.*(velocity_vec_length + a_qp))/hvel_qp, 2.) + pow((4.*_mu_qp)/(density*pow(hvel_qp, 2.)), 2.), -1./2.);      // NOTE: assumes dtime = 1. [Steady-State]
              
              // Calculate U-Momentum Navier-Stokes Solution for Stabilization
              NS_stab_u = (a1_urow.add(1., b1_urow)).dot(dUdx) + (a2_urow.add(1., b2_urow)).dot(dUdy) + (a3_urow.add(1., b3_urow)).dot(dUdz) 
                          - c11_urow.dot(momementum_gradphi[ii][qp](0)*dUdx)
                          - c12_urow.dot(momementum_gradphi[ii][qp](0)*dUdy)
                          - c13_urow.dot(momementum_gradphi[ii][qp](0)*dUdz)
                          - c21_urow.dot(momementum_gradphi[ii][qp](1)*dUdx)
                          - c22_urow.dot(momementum_gradphi[ii][qp](1)*dUdy)
                          - c23_urow.dot(momementum_gradphi[ii][qp](1)*dUdz)
                          - c31_urow.dot(momementum_gradphi[ii][qp](2)*dUdx)
                          - c32_urow.dot(momementum_gradphi[ii][qp](2)*dUdy)
                          - c33_urow.dot(momementum_gradphi[ii][qp](2)*dUdz);
              
              // F{rho_u}
              Frho_u(ii) -= JxW_momentum[qp] * 
                    // Conservative Navier-Stokes
                    (momentum_phi[ii][qp] * (a1_urow.dot(dUdx) + a2_urow.dot(dUdy) + a3_urow.dot(dUdz)) +
                    momentum_gradphi[ii][qp](0) * (c11_urow.dot(dUdx) + c12_urow.dot(dUdy) + c13_urow.dot(dUdz)) +
                    momentum_gradphi[ii][qp](1) * (c21_urow.dot(dUdx) + c22_urow.dot(dUdy) + c23_urow.dot(dUdz)) +
                    momentum_gradphi[ii][qp](2) * (c31_urow.dot(dUdx) + c32_urow.dot(dUdy) + c33_urow.dot(dUdz)) 
                    // SUPG Stabilization
                    + stab_SUPG_momentum * 
                    ((momentum_gradphi[ii][qp](0) * a1_urow.dot(NS_stab_u)) +
                     (momentum_gradphi[ii][qp](1) * a2_urow.dot(NS_stab_u)) +
                     (momentum_gradphi[ii][qp](2) * a3_urow.dot(NS_stab_u)))
                    // Shock Capturing Operator -- IN Process
                    );
                    
              if(std::isnan(Frho_u(ii)))
                { std::cout << "--- assemble_momentum_energy_time_derivative() ---" << "\n";
                  std::cout << "Frho_u(ii) = " << Frho_u(ii) << "\n";
                  std::cout << "ii = " << ii << "\n";
                  std::cout << "JxW_momentum[qp] = " << JxW_momentum[qp] << "\n";
                  std::cout << "a1_urow.dot(dUdx) = " << a1_urow.dot(dUdx) << "\n";
                  std::cout << "a2_urow.dot(dUdy) = " << a2_urow.dot(dUdy) << "\n";
                  std::cout << "a3_urow.dot(dUdz) = " << a3_urow.dot(dUdz) << "\n";
                  std::cout << "\n";
                  std::cout << "density = " << density << "\n";
                  std::cout << "\n";
                  std::cout << "dUdx(0) = " << dUdx(0) << "\n";
                  std::cout << "dUdx(1) = " << dUdx(1) << "\n"; 
                  std::cout << "dUdx(2) = " << dUdx(2) << "\n"; 
                  std::cout << "dUdx(3) = " << dUdx(3) << "\n"; 
                  std::cout << "dUdx(4) = " << dUdx(4) << "\n";
                  std::cout << "dUdy(0) = " << dUdy(0) << "\n";
                  std::cout << "dUdy(1) = " << dUdy(1) << "\n"; 
                  std::cout << "dUdy(2) = " << dUdy(2) << "\n"; 
                  std::cout << "dUdy(3) = " << dUdy(3) << "\n"; 
                  std::cout << "dUdy(4) = " << dUdy(4) << "\n";  
                  std::cout << "dUdz(0) = " << dUdz(0) << "\n";
                  std::cout << "dUdz(1) = " << dUdz(1) << "\n"; 
                  std::cout << "dUdz(2) = " << dUdz(2) << "\n"; 
                  std::cout << "dUdz(3) = " << dUdz(3) << "\n"; 
                  std::cout << "dUdz(4) = " << dUdz(4) << "\n";
                  
                  std::cout << "a1_urow(0) = " << a1_urow(0) << "\n";
                  std::cout << "a1_urow(1) = " << a1_urow(1) << "\n";
                  std::cout << "a1_urow(2) = " << a1_urow(2) << "\n";
                  std::cout << "a1_urow(3) = " << a1_urow(3) << "\n";
                  std::cout << "a1_urow(4) = " << a1_urow(4) << "\n";
                  std::cout << "a1_vrow(0) = " << a1_vrow(0) << "\n";
                  std::cout << "a1_vrow(1) = " << a1_vrow(1) << "\n";
                  std::cout << "a1_vrow(2) = " << a1_vrow(2) << "\n";
                  std::cout << "a1_vrow(3) = " << a1_vrow(3) << "\n";
                  std::cout << "a1_vrow(4) = " << a1_vrow(4) << "\n";
                  std::cout << "a1_wrow(0) = " << a1_wrow(0) << "\n";
                  std::cout << "a1_wrow(1) = " << a1_wrow(1) << "\n";
                  std::cout << "a1_wrow(2) = " << a1_wrow(2) << "\n";
                  std::cout << "a1_wrow(3) = " << a1_wrow(3) << "\n";
                  std::cout << "a1_wrow(4) = " << a1_wrow(4) << "\n";
                  
                  std::cout << "a2_urow(0) = " << a2_urow(0) << "\n";
                  std::cout << "a2_urow(1) = " << a2_urow(1) << "\n";
                  std::cout << "a2_urow(2) = " << a2_urow(2) << "\n";
                  std::cout << "a2_urow(3) = " << a2_urow(3) << "\n";
                  std::cout << "a2_urow(4) = " << a2_urow(4) << "\n";
                  std::cout << "a2_vrow(0) = " << a2_vrow(0) << "\n";
                  std::cout << "a2_vrow(1) = " << a2_vrow(1) << "\n";
                  std::cout << "a2_vrow(2) = " << a2_vrow(2) << "\n";
                  std::cout << "a2_vrow(3) = " << a2_vrow(3) << "\n";
                  std::cout << "a2_vrow(4) = " << a2_vrow(4) << "\n";
                  std::cout << "a2_wrow(0) = " << a2_wrow(0) << "\n";
                  std::cout << "a2_wrow(1) = " << a2_wrow(1) << "\n";
                  std::cout << "a2_wrow(2) = " << a2_wrow(2) << "\n";
                  std::cout << "a2_wrow(3) = " << a2_wrow(3) << "\n";
                  std::cout << "a2_wrow(4) = " << a2_wrow(4) << "\n";
                  
                  std::cout << "a3_urow(0) = " << a3_urow(0) << "\n";
                  std::cout << "a3_urow(1) = " << a3_urow(1) << "\n";
                  std::cout << "a3_urow(2) = " << a3_urow(2) << "\n";
                  std::cout << "a3_urow(3) = " << a3_urow(3) << "\n";
                  std::cout << "a3_urow(4) = " << a3_urow(4) << "\n";
                  std::cout << "a3_vrow(0) = " << a3_vrow(0) << "\n";
                  std::cout << "a3_vrow(1) = " << a3_vrow(1) << "\n";
                  std::cout << "a3_vrow(2) = " << a3_vrow(2) << "\n";
                  std::cout << "a3_vrow(3) = " << a3_vrow(3) << "\n";
                  std::cout << "a3_vrow(4) = " << a3_vrow(4) << "\n";
                  std::cout << "a3_wrow(0) = " << a3_wrow(0) << "\n";
                  std::cout << "a3_wrow(1) = " << a3_wrow(1) << "\n";
                  std::cout << "a3_wrow(2) = " << a3_wrow(2) << "\n";
                  std::cout << "a3_wrow(3) = " << a3_wrow(3) << "\n";
                  std::cout << "a3_wrow(4) = " << a3_wrow(4) << "\n";
                  
                  std::cout << " -----------------------" << "\n";
                }

              // Calculate V-Momentum Navier-Stokes Solution for Stabilization
              NS_stab_v = (a1_vrow.add(1., b1_vrow)).dot(dUdx) + (a2_vrow.add(1., b2_vrow)).dot(dUdy) + (a3_vrow.add(1., b3_vrow)).dot(dUdz) 
                          - c11_vrow.dot(momementum_gradphi[ii][qp](0)*dUdx)
                          - c12_vrow.dot(momementum_gradphi[ii][qp](0)*dUdy)
                          - c13_vrow.dot(momementum_gradphi[ii][qp](0)*dUdz)
                          - c21_vrow.dot(momementum_gradphi[ii][qp](1)*dUdx)
                          - c22_vrow.dot(momementum_gradphi[ii][qp](1)*dUdy)
                          - c23_vrow.dot(momementum_gradphi[ii][qp](1)*dUdz)
                          - c31_vrow.dot(momementum_gradphi[ii][qp](2)*dUdx)
                          - c32_vrow.dot(momementum_gradphi[ii][qp](2)*dUdy)
                          - c33_vrow.dot(momementum_gradphi[ii][qp](2)*dUdz);
              
              // F{rho_v}     
              Frho_v(ii) -= JxW_momentum[qp] * 
                    // Conservative Navier-Stokes
                    (momentum_phi[ii][qp] * (a1_vrow.dot(dUdx) + a2_vrow.dot(dUdy) + a3_vrow.dot(dUdz)) + 
                    momentum_gradphi[ii][qp](0) * (c11_vrow.dot(dUdx) + c12_vrow.dot(dUdy) + c13_vrow.dot(dUdz)) +
                    momentum_gradphi[ii][qp](1) * (c21_vrow.dot(dUdx) + c22_vrow.dot(dUdy) + c23_vrow.dot(dUdz)) +
                    momentum_gradphi[ii][qp](2) * (c31_vrow.dot(dUdx) + c32_vrow.dot(dUdy) + c33_vrow.dot(dUdz)) 
                    // SUPG Stabilization
                    + stab_SUPG_momentum * 
                    ((momentum_gradphi[ii][qp](0) * a1_vrow.dot(NS_stab_v)) +
                     (momentum_gradphi[ii][qp](1) * a2_vrow.dot(NS_stab_v)) +
                     (momentum_gradphi[ii][qp](2) * a3_vrow.dot(NS_stab_v)))
                    // Shock Capturing Operator -- IN Process
                    );
                    
              if(std::isnan(Frho_v(ii)))
                { std::cout << "--- assemble_momentum_energy_time_derivative() ---" << "\n";
                  std::cout << "Frho_v(ii) = " << Frho_v(ii) << "\n";
                  std::cout << "ii = " << ii << "\n";
                  std::cout << "JxW_momentum[qp] = " << JxW_momentum[qp] << "\n";
                  std::cout << "a1_vrow.dot(dUdx) = " << a1_vrow.dot(dUdx) << "\n";
                  std::cout << "a2_vrow.dot(dUdy) = " << a2_vrow.dot(dUdy) << "\n";
                  std::cout << "a3_vrow.dot(dUdz) = " << a3_vrow.dot(dUdz) << "\n";
                  std::cout << " -----------------------" << "\n";
                }
                    
              if (this->_momentum_vars.dim() == 3)
                {
                  // Calculate W-Momentum Navier-Stokes Solution for Stabilization
                  NS_stab_w = (a1_wrow.add(1., b1_wrow)).dot(dUdx) + (a2_wrow.add(1., b2_wrow)).dot(dUdy) + (a3_wrow.add(1., b3_wrow)).dot(dUdz) 
                          - c11_wrow.dot(momementum_gradphi[ii][qp](0)*dUdx)
                          - c12_wrow.dot(momementum_gradphi[ii][qp](0)*dUdy)
                          - c13_wrow.dot(momementum_gradphi[ii][qp](0)*dUdz)
                          - c21_wrow.dot(momementum_gradphi[ii][qp](1)*dUdx)
                          - c22_wrow.dot(momementum_gradphi[ii][qp](1)*dUdy)
                          - c23_wrow.dot(momementum_gradphi[ii][qp](1)*dUdz)
                          - c31_wrow.dot(momementum_gradphi[ii][qp](2)*dUdx)
                          - c32_wrow.dot(momementum_gradphi[ii][qp](2)*dUdy)
                          - c33_wrow.dot(momementum_gradphi[ii][qp](2)*dUdz);
                
                
                  // F{rho_w}
                  (*Frho_w)(ii) -= JxW_momentum[qp] *
                      // Conservative Navier-Stokes
                      (momentum_phi[ii][qp] * (a1_wrow.dot(dUdx) + a2_wrow.dot(dUdy) + a3_wrow.dot(dUdz)) + 
                      momentum_gradphi[ii][qp](0) * (c11_wrow.dot(dUdx) + c12_wrow.dot(dUdy) + c13_wrow.dot(dUdz)) +
                      momentum_gradphi[ii][qp](1) * (c21_wrow.dot(dUdx) + c22_wrow.dot(dUdy) + c23_wrow.dot(dUdz)) +
                      momentum_gradphi[ii][qp](2) * (c31_wrow.dot(dUdx) + c32_wrow.dot(dUdy) + c33_wrow.dot(dUdz))
                      // SUPG Stabilization
                      + stab_SUPG_momentum * 
                      ((momentum_gradphi[ii][qp](0) * a1_wrow.dot(NS_stab_w)) +
                       (momentum_gradphi[ii][qp](1) * a2_wrow.dot(NS_stab_w)) +
                       (momentum_gradphi[ii][qp](2) * a3_wrow.dot(NS_stab_w)))
                      // Shock Capturing Operator -- IN Process                      
                      );
                      
                  if(std::isnan((*Frho_w)(ii)))
                    { std::cout << "--- assemble_momentum_energy_time_derivative() ---" << "\n";
                      std::cout << "Frho_w(ii) = " << (*Frho_w)(ii) << "\n";
                      std::cout << "ii = " << ii << "\n";
                      std::cout << "JxW_momentum[qp] = " << JxW_momentum[qp] << "\n";
                      std::cout << " -----------------------" << "\n";
                    }                
                }  // End of 3D if statment 
                
              /* if (compute_jacobian)
                {
                  // Krho_u_rho, Krho_v_rho, Krho_w_rho
                  for (unsigned int jj=0; jj != n_rho_dofs; jj++)
                    {
                      // --- NOTE: May be able to reuse these vectors in all K(ii,jj) calculation to save on memory  --- //
                      d_dx_rho[0] = rho_gradphi[jj][qp](0);          // dphi/dx for rho all other d_dx_rho[] = 0. from initialization
                      d_dy_rho[0] = rho_gradphi[jj][qp](1);          // dphi/dy for rho all other d_dy_rho[] = 0. from initialization
                      d_dz_rho[0] = rho_gradphi[jj][qp](2);          // dphi/dz for rho all other d_dz_rho[] = 0. from initialization
                      
                      Krho_u_rho(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                          (momentum_phi[ii][qp] * ((a1_urow + b1_urow)*d_dx_rho + (a2_urow + b2_urow)*d_dy_rho + (a3_urow + b3_urow)*d_dz_rho) +    // inviscid & inviscid-viscid interaction portions
                          momentum_gradphi[ii][qp](0) * (c11_urow*d_dx_rho + c12_urow*d_dy_rho + c13_urow*d_dz_rho) +
                          momentum_gradphi[ii][qp](1) * (c21_urow*d_dx_rho + c22_urow*d_dy_rho + c23_urow*d_dz_rho) +
                          momentum_gradphi[ii][qp](2) * (c31_urow*d_dx_rho + c32_urow*d_dy_rho + c33_urow*d_dz_rho)
                          );
                          
                      Krho_v_rho(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                          (momentum_phi[ii][qp] * ((a1_vrow + b1_vrow)*d_dx_rho + (a2_vrow + b2_vrow)*d_dy_rho + (a3_vrow + b3_vrow)*d_dz_rho) +    // inviscid & inviscid-viscid interaction portions
                          momentum_gradphi[ii][qp](0) * (c11_vrow*d_dx_rho + c12_vrow*d_dy_rho + c13_vrow*d_dz_rho) +
                          momentum_gradphi[ii][qp](1) * (c21_vrow*d_dx_rho + c22_vrow*d_dy_rho + c23_vrow*d_dz_rho) +
                          momentum_gradphi[ii][qp](2) * (c31_vrow*d_dx_rho + c32_vrow*d_dy_rho + c33_vrow*d_dz_rho)
                          );
                          
                      if (this->_momentum_vars.dim() == 3)
                        {
                          (*Krho_w_rho)(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                              (momentum_phi[ii][qp] * ((a1_wrow + b1_wrow)*d_dx_rho + (a2_wrow + b2_wrow)*d_dy_rho + (a3_wrow + b3_wrow)*d_dz_rho) +    // inviscid & inviscid-viscid interaction portions
                              momentum_gradphi[ii][qp](0) * (c11_wrow*d_dx_rho + c12_wrow*d_dy_rho + c13_wrow*d_dz_rho) +
                              momentum_gradphi[ii][qp](1) * (c21_wrow*d_dx_rho + c22_wrow*d_dy_rho + c23_wrow*d_dz_rho) +
                              momentum_gradphi[ii][qp](2) * (c31_wrow*d_dx_rho + c32_wrow*d_dy_rho + c33_wrow*d_dz_rho)
                              );
                        }  // End of 3D Krho_v_rho if statement for 3D flows
                    }  // End of Krho_u_rho, Krho_v_rho Krho_w_rho Jacobian Loop
                    
                  // Calculate Krho_u_rho_u, Krho_u_rho_v, etc. Jacobians
                  for (unsigned int jj=0; jj != n_rho_u_dofs; jj++)
                    {
                      // --- NOTE: May be be able to find a way to reuse d_dx_ABC vectros to save on memory  --- //
                      d_dx_umomentum[1] = momentum_gradphi[jj][qp](0);
                      d_dy_umomentum[1] = momentum_gradphi[jj][qp](1);
                      d_dz_umomentum[1] = momentum_gradphi[jj][qp](2);
                      
                      d_dx_vmomentum[2] = d_dx_umomentum[1];
                      d_dy_vmomentum[2] = d_dy_umomentum[1];
                      d_dz_vmomentum[2] = d_dz_umomentum[1];
                      
                      if (this->_momentum_vars.dim() == 3)
                        {
                          d_dx_wmomentum[3] = d_dx_umomentum[1];
                          d_dy_wmomentum[3] = d_dy_umomentum[1];
                          d_dz_wmomentum[3] = d_dz_umomentum[1];
                        }
                      
                      Krho_u_rho_u(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                          (momentum_phi[ii][qp] * ((a1_urow + b1_urow)*d_dx_umomentum + (a2_urow + b2_urow)*d_dy_umomentum + (a3_urow + b3_urow)*d_dz_umomentum) +    // inviscid & inviscid-viscid interaction portions
                          momentum_gradphi[ii][qp](0) * (c11_urow*d_dx_umomentum + c12_urow*d_dy_umomentum + c13_urow*d_dz_umomentum) +
                          momentum_gradphi[ii][qp](1) * (c21_urow*d_dx_umomentum + c22_urow*d_dy_umomentum + c23_urow*d_dz_umomentum) +
                          momentum_gradphi[ii][qp](2) * (c31_urow*d_dx_umomentum + c32_urow*d_dy_umomentum + c33_urow*d_dz_umomentum)
                          );
                          
                      Krho_u_rho_v(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                          (momentum_phi[ii][qp] * ((a1_urow + b1_urow)*d_dx_vmomentum + (a2_urow + b2_urow)*d_dy_vmomentum + (a3_urow + b3_urow)*d_dz_vmomentum) +    // inviscid & inviscid-viscid interaction portions
                          momentum_gradphi[ii][qp](0) * (c11_urow*d_dx_vmomentum + c12_urow*d_dy_vmomentum + c13_urow*d_dz_vmomentum) +
                          momentum_gradphi[ii][qp](1) * (c21_urow*d_dx_vmomentum + c22_urow*d_dy_vmomentum + c23_urow*d_dz_vmomentum) +
                          momentum_gradphi[ii][qp](2) * (c31_urow*d_dx_vmomentum + c32_urow*d_dy_vmomentum + c33_urow*d_dz_vmomentum)
                          );
                          
                      Krho_v_rho_u(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                          (momentum_phi[ii][qp] * ((a1_vrow + b1_vrow)*d_dx_umomentum + (a2_vrow + b2_vrow)*d_dy_umomentum + (a3_vrow + b3_vrow)*d_dz_umomentum) +    // inviscid & inviscid-viscid interaction portions
                          momentum_gradphi[ii][qp](0) * (c11_vrow*d_dx_umomentum + c12_vrow*d_dy_umomentum + c13_vrow*d_dz_umomentum) +
                          momentum_gradphi[ii][qp](1) * (c21_vrow*d_dx_umomentum + c22_vrow*d_dy_umomentum + c23_vrow*d_dz_umomentum) +
                          momentum_gradphi[ii][qp](2) * (c31_vrow*d_dx_umomentum + c32_vrow*d_dy_umomentum + c33_vrow*d_dz_umomentum)
                          );
                      
                      Krho_v_rho_v(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                          (momentum_phi[ii][qp] * ((a1_vrow + b1_vrow)*d_dx_vmomentum + (a2_vrow + b2_vrow)*d_dy_vmomentum + (a3_vrow + b3_vrow)*d_dz_vmomentum) +    // inviscid & inviscid-viscid interaction portions
                          momentum_gradphi[ii][qp](0) * (c11_vrow*d_dx_vmomentum + c12_vrow*d_dy_vmomentum + c13_vrow*d_dz_vmomentum) +
                          momentum_gradphi[ii][qp](1) * (c21_vrow*d_dx_vmomentum + c22_vrow*d_dy_vmomentum + c23_vrow*d_dz_vmomentum) +
                          momentum_gradphi[ii][qp](2) * (c31_vrow*d_dx_vmomentum + c32_vrow*d_dy_vmomentum + c33_vrow*d_dz_vmomentum)
                          );
                          
                      if (this->_momentum_vars.dim() == 3)
                        {
                          (*Krho_u_rho_w)(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                              (momentum_phi[ii][qp] * ((a1_urow + b1_urow)*d_dx_wmomentum + (a2_urow + b2_urow)*d_dy_wmomentum + (a3_urow + b3_urow)*d_dz_wmomentum) +    // inviscid & inviscid-viscid interaction portions
                              momentum_gradphi[ii][qp](0) * (c11_urow*d_dx_wmomentum + c12_urow*d_dy_wmomentum + c13_urow*d_dz_wmomentum) +
                              momentum_gradphi[ii][qp](1) * (c21_urow*d_dx_wmomentum + c22_urow*d_dy_wmomentum + c23_urow*d_dz_wmomentum) +
                              momentum_gradphi[ii][qp](2) * (c31_urow*d_dx_wmomentum + c32_urow*d_dy_wmomentum + c33_urow*d_dz_wmomentum)
                              );
                              
                          (*Krho_v_rho_w)(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                              (momentum_phi[ii][qp] * ((a1_vrow + b1_vrow)*d_dx_wmomentum + (a2_vrow + b2_vrow)*d_dy_wmomentum + (a3_vrow + b3_vrow)*d_dz_wmomentum) +    // inviscid & inviscid-viscid interaction portions
                              momentum_gradphi[ii][qp](0) * (c11_vrow*d_dx_wmomentum + c12_vrow*d_dy_wmomentum + c13_vrow*d_dz_wmomentum) +
                              momentum_gradphi[ii][qp](1) * (c21_vrow*d_dx_wmomentum + c22_vrow*d_dy_wmomentum + c23_vrow*d_dz_wmomentum) +
                              momentum_gradphi[ii][qp](2) * (c31_vrow*d_dx_wmomentum + c32_vrow*d_dy_wmomentum + c33_vrow*d_dz_wmomentum)
                              );
                        
                          (*Krho_w_rho_u)(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                              (momentum_phi[ii][qp] * ((a1_wrow + b1_wrow)*d_dx_umomentum + (a2_wrow + b2_wrow)*d_dy_umomentum + (a3_wrow + b3_wrow)*d_dz_umomentum) +    // inviscid & inviscid-viscid interaction portions
                              momentum_gradphi[ii][qp](0) * (c11_wrow*d_dx_umomentum + c12_wrow*d_dy_umomentum + c13_wrow*d_dz_umomentum) +
                              momentum_gradphi[ii][qp](1) * (c21_wrow*d_dx_umomentum + c22_wrow*d_dy_umomentum + c23_wrow*d_dz_umomentum) +
                              momentum_gradphi[ii][qp](2) * (c31_wrow*d_dx_umomentum + c32_wrow*d_dy_umomentum + c33_wrow*d_dz_umomentum)
                              );
                              
                          (*Krho_w_rho_v)(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                              (momentum_phi[ii][qp] * ((a1_wrow + b1_wrow)*d_dx_vmomentum + (a2_wrow + b2_wrow)*d_dy_vmomentum + (a3_wrow + b3_wrow)*d_dz_vmomentum) +    // inviscid & inviscid-viscid interaction portions
                              momentum_gradphi[ii][qp](0) * (c11_wrow*d_dx_vmomentum + c12_wrow*d_dy_vmomentum + c13_wrow*d_dz_vmomentum) +
                              momentum_gradphi[ii][qp](1) * (c21_wrow*d_dx_vmomentum + c22_wrow*d_dy_vmomentum + c23_wrow*d_dz_vmomentum) +
                              momentum_gradphi[ii][qp](2) * (c31_wrow*d_dx_vmomentum + c32_wrow*d_dy_vmomentum + c33_wrow*d_dz_vmomentum)
                              );
                              
                          (*Krho_w_rho_w)(ii, jj) -= JxW * context.get_elem_solution_derivative() *
                              (momentum_phi[ii][qp] * ((a1_wrow + b1_wrow)*d_dx_wmomentum + (a2_wrow + b2_wrow)*d_dy_wmomentum + (a3_wrow + b3_wrow)*d_dz_wmomentum) +    // inviscid & inviscid-viscid interaction portions
                              momentum_gradphi[ii][qp](0) * (c11_wrow*d_dx_wmomentum + c12_wrow*d_dy_wmomentum + c13_wrow*d_dz_wmomentum) +
                              momentum_gradphi[ii][qp](1) * (c21_wrow*d_dx_wmomentum + c22_wrow*d_dy_wmomentum + c23_wrow*d_dz_wmomentum) +
                              momentum_gradphi[ii][qp](2) * (c31_wrow*d_dx_wmomentum + c32_wrow*d_dy_wmomentum + c33_wrow*d_dz_wmomentum)
                              );        
                        }  // End of 3D Krho_u_rho_w, etc if statement for 3D flows
                      
                    }  // End of Jacobian for loop
                
                }  // End of compute_jacobian IF statement
                */
            }  // End of ii loop for velocity degree of freedom
            
          for (unsigned int ii=0; ii != n_conserv_energy_dofs; ii++)
            {
              // Calculate Flow Aligned Length Scale
              hvel_qp = 2./(abs(unit_velocity.dot(conserv_energy_gradphi[ii][qp])));
              
              // Calculate SUPG Momentum Stabilization Factor
              stab_SUPG_energy = pow(pow(2./dtime, 2.) + pow((2.*(velocity_vec_length + a_qp))/hvel_qp, 2.) + pow((4.*_k_qp)/(density*_cp_qp*pow(hvel_qp, 2.)), 2.), -1./2.);      // NOTE: assumes dtime = 1. [Steady-State]
              
              // Calculate Conservative Energy Navier-Stokes Solution for Stabilization
              NS_stab_energy = (a1_energyrow.add(1., b1_energyrow)).dot(dUdx) + (a2_energyrow.add(1., b2_energyrow)).dot(dUdy) + (a3_energyrow.add(1., b3_energyrow)).dot(dUdz) 
                              - c11_energyrow.dot(conserv_energy_gradphi[ii][qp](0)*dUdx)
                              - c12_energyrow.dot(conserv_energy_gradphi[ii][qp](0)*dUdy)
                              - c13_energyrow.dot(conserv_energy_gradphi[ii][qp](0)*dUdz)
                              - c21_energyrow.dot(conserv_energy_gradphi[ii][qp](1)*dUdx)
                              - c22_energyrow.dot(conserv_energy_gradphi[ii][qp](1)*dUdy)
                              - c23_energyrow.dot(conserv_energy_gradphi[ii][qp](1)*dUdz)
                              - c31_energyrow.dot(conserv_energy_gradphi[ii][qp](2)*dUdx)
                              - c32_energyrow.dot(conserv_energy_gradphi[ii][qp](2)*dUdy)
                              - c33_energyrow.dot(conserv_energy_gradphi[ii][qp](2)*dUdz);
            
              // F{conserv_energy}
              Fconserv_energy(ii) -= JxW_energy[qp] *
                    // Conservative Navier-Stokes
                    (conserv_energy_phi[ii][qp] * (a1_energyrow.dot(dUdx) + a2_energyrow.dot(dUdy) + a3_energyrow.dot(dUdz)) +
                    conserv_energy_gradphi[ii][qp](0) * (c11_energyrow.dot(dUdx) + c12_energyrow.dot(dUdy) + c13_energyrow.dot(dUdz)) +
                    conserv_energy_gradphi[ii][qp](1) * (c21_energyrow.dot(dUdx) + c22_energyrow.dot(dUdy) + c23_energyrow.dot(dUdz)) +
                    conserv_energy_gradphi[ii][qp](2) * (c31_energyrow.dot(dUdx) + c32_energyrow.dot(dUdy) + c33_energyrow.dot(dUdz))
                    // SUPG Stabilization
                    + stab_SUPG_energy * 
                    ((conserv_energy_gradphi[ii][qp](0) * a1_energyrow.dot(NS_stab_energy)) +
                     (conserv_energy_gradphi[ii][qp](1) * a2_energyrow.dot(NS_stab_energy)) +
                     (conserv_energy_gradphi[ii][qp](2) * a3_energyrow.dot(NS_stab_energy)))
                    // Shock Capturing Operator -- IN Process  
                    );    
                    
              if(std::isnan(Fconserv_energy(ii)))
                { std::cout << "--- assemble_momentum_energy_time_derivative() ---" << "\n";
                  std::cout << "Fconserv_energy(ii) = " << Fconserv_energy(ii) << "\n";
                  std::cout << "ii = " << ii << "\n";
                  std::cout << "JxW_energy[qp] = " << JxW_energy[qp] << "\n";
                  std::cout << " -----------------------" << "\n";
                }        
            }  // End of ii for energy degree of freedom
        }    // End of Quadrature Point For Loop
      
      return;            
    }    // End of assemble_momentum_time_derivative()
    
    // --- Register Parameters
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::register_parameter( const std::string & param_name,
                                                                   libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer)
                                                                   const
    {
      ParameterUser::register_parameter( param_name, param_pointer);
      _mu.register_parameter(param_name, param_pointer);
      _cp.register_parameter(param_name, param_pointer);
      _k.register_parameter(param_name, param_pointer);
    }
    
    // --- Compute Postprocessed Quantities [MAY NOT NEED THIS? MAYBE WHERE I CALCULATE OUTPUT VALUES FOR SOLUTION??]
    template<class Mu, class SH, class TC>
    void ConservativeNavierStokes<Mu, SH, TC>::compute_postprocessed_quantity( unsigned int quantity_index,
                                                                               const AssemblyContext& context,
                                                                               const libMesh::Point& point,
                                                                               libMesh::Real& value )
    {
      // ----------------------------------------------------
      // ----- Add to code to calculate values of       -----
      // ----- variables required for the next          -----
      // ----- iteration?/step?                         -----
      // ----- May need this for Pressure/Temperature   -----
      // ----- in Compressible Navier Stokes Solutions  -----
      // ----------------------------------------------------
    }
         
}  // namespace GRINS

// Instantiate
template class GRINS::ConservativeNavierStokes<GRINS::ConstantViscosity,GRINS::ConstantSpecificHeat,GRINS::ConstantConductivity>;