! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   libAtoms+QUIP: atomistic simulation library
! HND X
! HND X   Portions of this code were written by
! HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! HND X
! HND X   Copyright 2006-2010.
! HND X
! HND X   Not for distribution
! HND X
! HND X   Portions of this code were written by Noam Bernstein as part of
! HND X   his employment for the U.S. Government, and are not subject
! HND X   to copyright in the USA.
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   http://www.libatoms.org
! HND X
! HND X  Additional contributions by
! HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module descriptors_module

   use error_module
   use system_module, only : dp, print, optional_default, system_timer, operator(//), split_string, string_to_int, split_string_simple, inoutput, OUTPUT, PRINT_VERBOSE, PRINT_NERD
   use linkedlist_module
   use units_module
   use linearalgebra_module
   use dictionary_module
   use paramreader_module
   use atoms_module
   use atoms_types_module
   use topology_module
   use mpi_context_module
   use table_module
   use permutation_maker_module
   use CInOutput_module
   use clusters_module
   use connection_module
   use angular_functions_module

   implicit none

   private
#ifdef GAP_VERSION
   integer, parameter :: gap_version = GAP_VERSION
#else
   integer, parameter :: gap_version = 0
#endif


   integer, parameter, public :: DT_NONE            =  0
   integer, parameter, public :: DT_BISPECTRUM_SO4  =  1
   integer, parameter, public :: DT_BISPECTRUM_SO3  =  2
   integer, parameter, public :: DT_BEHLER          =  3
   integer, parameter, public :: DT_DISTANCE_2B     =  4
   integer, parameter, public :: DT_COORDINATION    =  5
   integer, parameter, public :: DT_ANGLE_3B        =  6
   integer, parameter, public :: DT_CO_ANGLE_3B     =  7
   integer, parameter, public :: DT_CO_DISTANCE_2B  =  8
   integer, parameter, public :: DT_COSNX           =  9
   integer, parameter, public :: DT_TRIHIS          = 10
   integer, parameter, public :: DT_WATER_MONOMER   = 11
   integer, parameter, public :: DT_WATER_DIMER     = 12
   integer, parameter, public :: DT_A2_DIMER        = 13
   integer, parameter, public :: DT_AB_DIMER        = 14
   integer, parameter, public :: DT_BOND_REAL_SPACE = 15
   integer, parameter, public :: DT_ATOM_REAL_SPACE = 16
   integer, parameter, public :: DT_POWER_SO3       = 17
   integer, parameter, public :: DT_POWER_SO4       = 18
   integer, parameter, public :: DT_SOAP            = 19
   integer, parameter, public :: DT_AN_MONOMER      = 20
   integer, parameter, public :: DT_GENERAL_MONOMER = 21
   integer, parameter, public :: DT_GENERAL_DIMER   = 22
   integer, parameter, public :: DT_GENERAL_TRIMER  = 23
   integer, parameter, public :: DT_RDF             = 24
   integer, parameter, public :: DT_AS_DISTANCE_2B  = 25
   integer, parameter, public :: DT_MOLECULE_LO_D   = 26
   integer, parameter, public :: DT_alex            = 27
   integer, parameter, public :: DT_COM_DIMER       = 28
   integer, parameter, public :: DT_DISTANCE_NB     = 29

   integer, parameter :: NP_WATER_DIMER    = 8
   integer, parameter :: NP_A2_DIMER       = 8
   integer, parameter :: NP_AB_DIMER       = 2

   type transfer_parameters_type
      logical :: do_transfer
      real(dp) :: factor, r0, width
   endtype transfer_parameters_type

   type descriptor_data_mono
      real(dp), dimension(:), allocatable :: data
      real(dp), dimension(:,:,:), allocatable :: grad_data
      ! ci : atom indices amongst which to distribute energy of descriptor
      ! ii : all atoms involved in descriptor (for partial derivatives)
      integer, dimension(:), allocatable :: ci, ii
      real(dp), dimension(:,:), allocatable :: pos
      logical :: has_data
      logical, dimension(:), allocatable :: has_grad_data

      real(dp) :: covariance_cutoff = 1.0_dp
      real(dp), dimension(:,:), allocatable :: grad_covariance_cutoff
   endtype descriptor_data_mono

   type cplx_2d
      complex(dp), dimension(:,:), allocatable :: mm
   endtype cplx_2d

   type real_2d
      real(dp), dimension(:,:), allocatable :: mm
   endtype real_2d

   type cplx_3d
      complex(dp), dimension(:,:,:), allocatable :: mm
   endtype cplx_3d

   type RadialFunction_type
      integer :: n_max
      real(dp) :: cutoff, min_cutoff
      real(dp), dimension(:,:), allocatable :: RadialTransform
      real(dp), dimension(:), allocatable :: NormFunction

      logical :: initialised = .false.
   endtype RadialFunction_type

   type fourier_SO4_type
      real(dp) :: cutoff
      real(dp) :: z0_ratio
      real(dp) :: z0
      integer :: j_max, Z
      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.
   endtype fourier_SO4_type

   type bispectrum_SO4
      real(dp), pointer :: cutoff
      integer, pointer :: j_max, Z
      real(dp), pointer :: z0_ratio
      real(dp), pointer :: z0

      integer, dimension(:), pointer :: species_Z
      real(dp), dimension(:), pointer :: w
      
      type(fourier_SO4_type) :: fourier_SO4

      logical :: initialised = .false.

   endtype bispectrum_SO4

   type bispectrum_SO3

      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.

   endtype bispectrum_SO3

   type behler_g2
      real(dp) :: eta
      real(dp) :: rs
      real(dp) :: rc
   endtype behler_g2

   type behler_g3
      real(dp) :: eta
      real(dp) :: lambda
      real(dp) :: zeta
      real(dp) :: rc
   endtype behler_g3

   type behler

      real(dp) :: cutoff = 0.0_dp
      logical :: initialised = .false.

      integer :: n_g2, n_g3
      type(behler_g2), dimension(:), allocatable :: g2
      type(behler_g3), dimension(:), allocatable :: g3

   endtype behler

   type distance_2b
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: Z1, Z2
      character(STRING_LENGTH) :: resid_name
      logical :: only_intra, only_inter

      logical :: initialised = .false.

   endtype distance_2b

   type coordination
      real(dp) :: cutoff
      real(dp) :: transition_width
      integer :: Z

      logical :: initialised = .false.

   endtype coordination

   type angle_3b
      real(dp) :: cutoff
      integer :: Z, Z1, Z2

      logical :: initialised = .false.

   endtype angle_3b

   type co_angle_3b
      real(dp) :: cutoff
      real(dp) :: coordination_cutoff
      real(dp) :: coordination_transition_width
      integer :: Z, Z1, Z2

      logical :: initialised = .false.

   endtype co_angle_3b

   type co_distance_2b
      real(dp) :: cutoff
      real(dp) :: transition_width
      real(dp) :: coordination_cutoff
      real(dp) :: coordination_transition_width
      integer :: Z1, Z2

      logical :: initialised = .false.

   endtype co_distance_2b

   type cosnx

      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.

   endtype cosnx

   type trihis
      real(dp) :: cutoff
      integer :: n_gauss

      real(dp), dimension(:,:), allocatable :: gauss_centre
      real(dp), dimension(:,:), allocatable :: gauss_width

      logical :: initialised = .false.

   endtype trihis

   type water_monomer
      real(dp) :: cutoff

      logical :: initialised = .false.

   endtype water_monomer

   type water_dimer
      real(dp) :: cutoff, cutoff_transition_width
      real(dp) :: monomer_cutoff
      logical :: OHH_ordercheck

      logical :: initialised = .false.

   endtype water_dimer

   type A2_dimer
      real(dp) :: cutoff
      real(dp) :: monomer_cutoff
      integer :: atomic_number

      logical :: initialised = .false.

   endtype A2_dimer

   type AB_dimer
      real(dp) :: cutoff
      real(dp) :: monomer_cutoff
      integer :: atomic_number1, atomic_number2

      logical :: initialised = .false.

   endtype AB_dimer

   type bond_real_space
      real(dp) :: bond_cutoff
      real(dp) :: bond_transition_width
      real(dp) :: cutoff
      real(dp) :: transition_width
      real(dp) :: atom_sigma
      integer :: max_neighbours

      logical :: initialised = .false.

   endtype bond_real_space

   type atom_real_space
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: l_max
      real(dp) :: alpha
      real(dp) :: zeta

      logical :: initialised = .false.

   endtype atom_real_space

   type power_so3
      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.
   endtype power_so3

   type power_SO4
      real(dp), pointer :: cutoff
      integer, pointer :: j_max, Z
      real(dp), pointer :: z0_ratio
      real(dp), pointer :: z0

      integer, dimension(:), pointer :: species_Z
      real(dp), dimension(:), pointer :: w
      
      type(fourier_SO4_type) :: fourier_SO4

      logical :: initialised = .false.

   endtype power_SO4

   type soap
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      real(dp) :: alpha, atom_sigma, covariance_sigma0, central_weight

      integer :: l_max, n_max, n_Z, n_species
      integer, dimension(:), allocatable :: species_Z, Z
      real(dp), dimension(:), allocatable :: r_basis
      real(dp), dimension(:,:), allocatable :: transform_basis,cholesky_overlap_basis

      logical :: global = .false.
      logical :: central_reference_all_species = .false.
      logical :: diagonal_radial = .false.
      logical :: normalise = .true.
      logical :: initialised = .false.
   endtype soap
   public :: soap
 
   type AN_monomer
      real(dp) :: cutoff
      integer :: atomic_number
      integer :: N

      logical :: initialised = .false.
      logical :: do_atomic = .false.

   endtype AN_monomer

   type general_monomer
      type(permutation_data_type) :: permutation_data
      integer, dimension(:), allocatable :: signature
      real(dp) :: cutoff, cutoff_transition_width
      logical :: atom_ordercheck, internal_swaps_only
      logical :: strict
      logical :: initialised = .false.
   endtype general_monomer
   public :: general_monomer

   type general_dimer
      type(permutation_data_type) :: permutation_data
      integer, dimension(:), allocatable :: signature_one, signature_two
      integer, dimension(:,:), allocatable :: component_atoms
      real(dp) :: cutoff, cutoff_transition_width, monomer_one_cutoff, monomer_two_cutoff
      logical :: atom_ordercheck, internal_swaps_only, use_smooth_cutoff, monomers_identical,double_count
      logical :: strict, use_com, mpifind
      type(transfer_parameters_type) :: transfer_parameters
      logical :: initialised = .false.
      logical, dimension(:), allocatable :: is_intermolecular, cutoff_contributor
   endtype general_dimer

   type general_trimer
      type(permutation_data_type) :: permutation_data
      integer, dimension(:), allocatable :: signature_one, signature_two, signature_three
      integer, dimension(:,:), allocatable :: component_atoms
      real(dp) :: cutoff, cutoff_transition_width, monomer_one_cutoff, monomer_two_cutoff, monomer_three_cutoff
      logical :: atom_ordercheck, internal_swaps_only, use_smooth_cutoff, one_two_identical, one_three_identical, two_three_identical
      logical :: strict, use_com, mpifind
      logical :: initialised = .false.
      logical, dimension(:), allocatable :: is_intermolecular, cutoff_contributor
   endtype general_trimer

   type rdf
      real(dp) :: cutoff
      real(dp) :: transition_width, w_gauss
      integer :: Z, n_gauss
      real(dp), dimension(:), allocatable :: r_gauss

      logical :: initialised = .false.

   endtype rdf

   type as_distance_2b
      real(dp) :: min_cutoff, max_cutoff, as_cutoff, overlap_alpha
      real(dp) :: min_transition_width, max_transition_width, as_transition_width
      real(dp) :: coordination_cutoff
      real(dp) :: coordination_transition_width
      integer :: Z1, Z2

      logical :: initialised = .false.

   endtype as_distance_2b

   type molecule_lo_d
      type(permutation_data_type) :: permutation_data
      type(Atoms) :: template_atoms
      integer :: n_atoms, max_dimension ! max_dimension is descriptor dimension if include all interatomic distances
      integer, dimension(:), allocatable :: signature, included_components
      integer, dimension(:,:), allocatable :: component_atoms
      real(dp) :: cutoff, cutoff_transition_width
      integer :: neighbour_graph_depth
      logical :: atom_ordercheck, use_smooth_cutoff
      logical :: initialised = .false.
      type(Table) :: bonds, atom_pairs
      integer :: desctype
   endtype molecule_lo_d

   type alex

      integer :: Z, power_min, power_max 
      real(dp) :: cutoff

      integer :: n_species
      integer, dimension(:), allocatable :: species_Z

      logical :: initialised = .false.
   endtype alex

   type com_dimer
      integer, dimension(:), allocatable :: signature_one, signature_two
      real(dp) :: cutoff, cutoff_transition_width, monomer_one_cutoff, monomer_two_cutoff
      logical :: atom_ordercheck, use_smooth_cutoff, monomers_identical
      logical :: strict, mpifind
      type(transfer_parameters_type) :: transfer_parameters
      logical :: initialised = .false.
      logical, dimension(:), allocatable :: is_intermolecular, cutoff_contributor
   endtype com_dimer

   type distance_Nb
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: order
      integer, dimension(:), allocatable :: Z
      integer :: n_permutations
      integer, dimension(:,:), allocatable :: permutations
      logical, dimension(:,:,:), allocatable :: monomerConnectivities
      logical :: compact_clusters = .false.
      logical :: initialised = .false.
   endtype distance_Nb

   type descriptor
      integer :: descriptor_type = DT_NONE

      type(bispectrum_SO4)  :: descriptor_bispectrum_SO4
      type(bispectrum_SO3)  :: descriptor_bispectrum_SO3
      type(behler)          :: descriptor_behler
      type(distance_2b)     :: descriptor_distance_2b
      type(coordination)    :: descriptor_coordination
      type(angle_3b)        :: descriptor_angle_3b
      type(co_angle_3b)     :: descriptor_co_angle_3b
      type(co_distance_2b)  :: descriptor_co_distance_2b
      type(cosnx)           :: descriptor_cosnx
      type(trihis)          :: descriptor_trihis
      type(water_monomer)   :: descriptor_water_monomer
      type(water_dimer)     :: descriptor_water_dimer
      type(A2_dimer)        :: descriptor_A2_dimer
      type(AB_dimer)        :: descriptor_AB_dimer
      type(bond_real_space) :: descriptor_bond_real_space
      type(atom_real_space) :: descriptor_atom_real_space
      type(power_so3)       :: descriptor_power_so3
      type(power_SO4)       :: descriptor_power_SO4
      type(soap)            :: descriptor_soap
      type(AN_monomer)      :: descriptor_AN_monomer
      type(general_monomer) :: descriptor_general_monomer
      type(general_dimer)   :: descriptor_general_dimer
      type(general_trimer)  :: descriptor_general_trimer
      type(rdf)             :: descriptor_rdf
      type(as_distance_2b)  :: descriptor_as_distance_2b
      type(molecule_lo_d)   :: descriptor_molecule_lo_d
      type(alex)           :: descriptor_alex
      type(com_dimer)       :: descriptor_com_dimer
      type(distance_Nb)     :: descriptor_distance_Nb
   endtype
   
   type descriptor_data
      type(descriptor_data_mono), dimension(:), allocatable :: x
   endtype descriptor_data
   
   type cplx_1d
      complex(dp), dimension(:), allocatable :: m
   endtype cplx_1d
   
   type real_1d
      real(dp), dimension(:), allocatable :: m
   endtype real_1d
   
   type spherical_harmonics_type
      type(cplx_1d), dimension(:), allocatable :: spherical_harmonics
      type(cplx_2d), dimension(:), allocatable :: grad_spherical_harmonics
      real(dp) :: r
      real(dp), dimension(3) :: u
   endtype spherical_harmonics_type
   
   type neighbour_type
      type(spherical_harmonics_type), dimension(:), allocatable :: neighbour
   endtype neighbour_type
   
   type grad_spherical_harmonics_overlap_type
      type(cplx_3d), dimension(:), allocatable :: grad_integral
   endtype grad_spherical_harmonics_overlap_type
   
   public :: neighbour_type, real_space_fourier_coefficients, real_space_covariance_coefficient
   public :: SphericalYCartesian
   
   interface initialise
      module procedure descriptor_initialise, RadialFunction_initialise, fourier_so4_initialise, &
      bispectrum_SO4_initialise, bispectrum_SO3_initialise, behler_initialise, distance_2b_initialise, &
      coordination_initialise, angle_3b_initialise, co_angle_3b_initialise, co_distance_2b_initialise, cosnx_initialise, trihis_initialise, &
      water_monomer_initialise, water_dimer_initialise, A2_dimer_initialise, AB_dimer_initialise, &
      bond_real_space_initialise, atom_real_space_initialise, power_so3_initialise, power_SO4_initialise, soap_initialise, AN_monomer_initialise, &
      general_monomer_initialise, general_dimer_initialise, general_trimer_initialise, rdf_initialise, as_distance_2b_initialise, molecule_lo_d_initialise, alex_initialise, &
      transfer_initialise, com_dimer_initialise, distance_Nb_initialise
   endinterface initialise
   public :: initialise

   interface finalise
      module procedure descriptor_finalise, descriptor_data_finalise, RadialFunction_finalise, fourier_so4_finalise, cplx_2d_array1_finalise, cplx_3d_array2_finalise, &
      bispectrum_SO4_finalise, bispectrum_SO3_finalise, behler_finalise, distance_2b_finalise, coordination_finalise, angle_3b_finalise, co_angle_3b_finalise, &
      co_distance_2b_finalise, cosnx_finalise, trihis_finalise, water_monomer_finalise, water_dimer_finalise, &
      A2_dimer_finalise, AB_dimer_finalise, bond_real_space_finalise, atom_real_space_finalise, power_so3_finalise, power_SO4_finalise, soap_finalise, &
      AN_monomer_finalise, general_monomer_finalise, general_dimer_finalise, general_trimer_finalise, rdf_finalise, as_distance_2b_finalise, molecule_lo_d_finalise, alex_finalise, com_dimer_finalise, &
      distance_Nb_finalise
   endinterface finalise
   public :: finalise

   interface calc
      module procedure descriptor_calc, descriptor_calc_array, bispectrum_SO4_calc, bispectrum_SO3_calc, behler_calc, distance_2b_calc, coordination_calc, angle_3b_calc, co_angle_3b_calc, &
      co_distance_2b_calc, cosnx_calc, trihis_calc, water_monomer_calc, water_dimer_calc, A2_dimer_calc, AB_dimer_calc, bond_real_space_calc, atom_real_space_calc, &
      power_so3_calc, power_SO4_calc, soap_calc, AN_monomer_calc, general_monomer_calc, general_dimer_calc, general_trimer_calc, rdf_calc, as_distance_2b_calc, molecule_lo_d_calc, alex_calc, com_dimer_calc, &
      distance_Nb_calc
   endinterface calc
   public :: calc

   interface cutoff
      module procedure descriptor_cutoff, bispectrum_SO4_cutoff, bispectrum_SO3_cutoff, behler_cutoff, distance_2b_cutoff, coordination_cutoff, angle_3b_cutoff, co_angle_3b_cutoff, &
      co_distance_2b_cutoff, cosnx_cutoff, trihis_cutoff, water_monomer_cutoff, water_dimer_cutoff, A2_dimer_cutoff, AB_dimer_cutoff, bond_real_space_cutoff, atom_real_space_cutoff, &
      power_so3_cutoff, power_SO4_cutoff, soap_cutoff, AN_monomer_cutoff, general_monomer_cutoff, general_dimer_cutoff, general_trimer_cutoff, rdf_cutoff, as_distance_2b_cutoff, &
      molecule_lo_d_cutoff, alex_cutoff, com_dimer_cutoff, distance_Nb_cutoff
   endinterface cutoff
   public :: cutoff

   interface descriptor_sizes
      module procedure descriptor_sizes, bispectrum_SO4_sizes, bispectrum_SO3_sizes, behler_sizes, distance_2b_sizes, coordination_sizes, angle_3b_sizes, co_angle_3b_sizes, &
      co_distance_2b_sizes, cosnx_sizes, trihis_sizes, water_monomer_sizes, water_dimer_sizes, A2_dimer_sizes, AB_dimer_sizes, bond_real_space_sizes, atom_real_space_sizes, &
      power_so3_sizes, power_SO4_sizes, soap_sizes, AN_monomer_sizes, general_monomer_sizes, general_dimer_sizes, general_trimer_sizes, rdf_sizes, as_distance_2b_sizes, &
      molecule_lo_d_sizes, alex_sizes, com_dimer_sizes, distance_Nb_sizes
   endinterface descriptor_sizes
   public :: descriptor_sizes

   public :: descriptor_MPI_setup

   public :: descriptor, descriptor_data, descriptor_dimensions, descriptor_n_permutations, descriptor_permutations, descriptor_str_add_species
   public :: real_space_covariance
   public :: cplx_1d, cplx_2d

   contains

   function get_descriptor_type(args_str,error)
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      integer :: get_descriptor_type

      type(Dictionary) :: params
      logical :: is_bispectrum_so4, is_bispectrum_so3, is_behler, is_distance_2b, is_coordination, is_angle_3b, &
         is_co_angle_3b, is_co_distance_2b, is_cosnx, is_trihis, is_water_monomer, is_water_dimer, is_A2_dimer, &
         is_AB_dimer, is_bond_real_space, is_atom_real_space, is_power_so3, is_power_so4, is_soap, &
         is_AN_monomer, is_general_monomer, is_general_dimer, is_general_trimer, is_rdf, is_as_distance_2b, &
         is_molecule_lo_d, is_alex, is_com_dimer, is_distance_Nb

      INIT_ERROR(error)

      call initialise(params)
      call param_register(params, 'bispectrum_so4', 'false', is_bispectrum_so4, help_string="Type of descriptor is bispectrum_so4.")
      call param_register(params, 'bispectrum_so3', 'false', is_bispectrum_so3, help_string="Type of descriptor is bispectrum_so3.")
      call param_register(params, 'behler', 'false', is_behler, help_string="Type of descriptor is behler.")
      call param_register(params, 'distance_2b', 'false', is_distance_2b, help_string="Type of descriptor is distance_2b.")
      call param_register(params, 'coordination', 'false', is_coordination, help_string="Type of descriptor is coordination.")
      call param_register(params, 'angle_3b', 'false', is_angle_3b, help_string="Type of descriptor is angle_3b.")
      call param_register(params, 'co_angle_3b', 'false', is_co_angle_3b, help_string="Type of descriptor is co_angle_3b.")
      call param_register(params, 'co_distance_2b', 'false', is_co_distance_2b, help_string="Type of descriptor is co_distance_2b.")
      call param_register(params, 'cosnx', 'false', is_cosnx, help_string="Type of descriptor is cosnx.")
      call param_register(params, 'trihis', 'false', is_trihis, help_string="Type of descriptor is trihis.")
      call param_register(params, 'water_monomer', 'false', is_water_monomer, help_string="Type of descriptor is water_monomer.")
      call param_register(params, 'water_dimer', 'false', is_water_dimer, help_string="Type of descriptor is water_dimer.")
      call param_register(params, 'A2_dimer', 'false', is_A2_dimer, help_string="Type of descriptor is A2_dimer.")
      call param_register(params, 'AB_dimer', 'false', is_AB_dimer, help_string="Type of descriptor is AB_dimer.")
      call param_register(params, 'bond_real_space', 'false', is_bond_real_space, help_string="Type of descriptor is bond_real_space.")
      call param_register(params, 'atom_real_space', 'false', is_atom_real_space, help_string="Type of descriptor is atom_real_space.")
      call param_register(params, 'power_so3', 'false', is_power_so3, help_string="Type of descriptor is power_so3.")
      call param_register(params, 'power_so4', 'false', is_power_so4, help_string="Type of descriptor is power_so4.")
      call param_register(params, 'soap', 'false', is_soap, help_string="Type of descriptor is soap.")
      call param_register(params, 'AN_monomer', 'false', is_AN_monomer, help_string="Type of descriptor is AN_monomer.")
      call param_register(params, 'general_monomer', 'false', is_general_monomer, help_string="Type of descriptor is general_monomer.")
      call param_register(params, 'general_dimer', 'false', is_general_dimer, help_string="Type of descriptor is general_dimer.")
      call param_register(params, 'general_trimer', 'false', is_general_trimer, help_string="Type of descriptor is general_trimer.")
      call param_register(params, 'rdf', 'false', is_rdf, help_string="Type of descriptor is rdf.")
      call param_register(params, 'as_distance_2b', 'false', is_as_distance_2b, help_string="Type of descriptor is as_distance_2b.")
      call param_register(params, 'molecule_lo_d', 'false', is_molecule_lo_d, help_string="Type of descriptor is molecule_lo_d.")
      call param_register(params, 'alex', 'false', is_alex, help_string="Type of descriptor is alex.")
      call param_register(params, 'com_dimer', 'false', is_com_dimer, help_string="Type of descriptor is com_dimer.")
      call param_register(params, 'distance_Nb', 'false', is_distance_Nb, help_string="Type of descriptor is distance_Nb.")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='descriptor_initialise args_str')) then
         RAISE_ERROR("descriptor_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if (count( (/is_bispectrum_so4, is_bispectrum_so3, is_behler, is_distance_2b, is_coordination, is_angle_3b, is_co_angle_3b, is_co_distance_2b, &
      is_cosnx, is_trihis, is_water_monomer, is_water_dimer, is_A2_dimer, is_AB_dimer, is_bond_real_space, is_atom_real_space, is_power_so3, is_power_so4, &
      is_soap, is_AN_monomer, is_general_monomer, is_general_dimer, is_general_trimer, is_rdf, is_as_distance_2b, is_molecule_lo_d, is_alex, is_com_dimer, &
      is_distance_Nb /) ) /= 1) then
         RAISE_ERROR("descriptor_initialise found too few or too many IP Model types args_str='"//trim(args_str)//"'", error)
      endif

      get_descriptor_type = DT_NONE

      if( is_bispectrum_so4 ) then
         get_descriptor_type = DT_BISPECTRUM_SO4
      elseif( is_bispectrum_so3 ) then
         get_descriptor_type = DT_BISPECTRUM_SO3
      elseif( is_behler ) then
         get_descriptor_type = DT_BEHLER
      elseif( is_distance_2b ) then
         get_descriptor_type = DT_DISTANCE_2B
      elseif( is_coordination ) then
         get_descriptor_type = DT_COORDINATION
      elseif( is_angle_3b ) then
         get_descriptor_type = DT_ANGLE_3B
      elseif( is_co_angle_3b ) then
         get_descriptor_type = DT_CO_ANGLE_3B
      elseif( is_co_distance_2b ) then
         get_descriptor_type = DT_CO_DISTANCE_2B
      elseif( is_cosnx ) then
         get_descriptor_type = DT_COSNX
      elseif( is_trihis ) then
         get_descriptor_type = DT_TRIHIS
      elseif( is_water_monomer ) then
         get_descriptor_type = DT_WATER_MONOMER
      elseif( is_water_dimer ) then
         get_descriptor_type = DT_WATER_DIMER
      elseif( is_A2_dimer ) then
         get_descriptor_type = DT_A2_DIMER
      elseif( is_AB_dimer ) then
         get_descriptor_type = DT_AB_DIMER
      elseif( is_bond_real_space ) then
         get_descriptor_type = DT_BOND_REAL_SPACE
      elseif( is_atom_real_space ) then
         get_descriptor_type = DT_ATOM_REAL_SPACE
      elseif( is_power_so3 ) then
         get_descriptor_type = DT_POWER_SO3
      elseif( is_power_so4 ) then
         get_descriptor_type = DT_POWER_SO4
      elseif( is_soap ) then
         get_descriptor_type = DT_SOAP
      elseif( is_AN_monomer ) then
         get_descriptor_type = DT_AN_MONOMER
      elseif( is_general_monomer ) then
         get_descriptor_type = DT_GENERAL_MONOMER
      elseif( is_general_dimer ) then
         get_descriptor_type = DT_GENERAL_DIMER
      elseif( is_general_trimer ) then
         get_descriptor_type = DT_GENERAL_TRIMER
      elseif( is_rdf ) then
         get_descriptor_type = DT_RDF
      elseif( is_as_distance_2b ) then
         get_descriptor_type = DT_AS_DISTANCE_2B
      elseif( is_molecule_lo_d ) then
         get_descriptor_type = DT_MOLECULE_LO_D
      elseif( is_alex ) then
         get_descriptor_type = DT_ALEX
      elseif( is_com_dimer ) then
         get_descriptor_type = DT_COM_DIMER
      elseif( is_distance_Nb ) then
         get_descriptor_type = DT_DISTANCE_NB
      endif

   endfunction get_descriptor_type

   subroutine descriptor_initialise(this,args_str,error)
      type(descriptor), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      this%descriptor_type = get_descriptor_type(args_str,error)

      select case(this%descriptor_type)
      case(DT_BISPECTRUM_SO4)
         call initialise(this%descriptor_bispectrum_SO4,args_str,error)
      case(DT_BISPECTRUM_SO3)
         call initialise(this%descriptor_bispectrum_SO3,args_str,error)
      case(DT_BEHLER)
         call initialise(this%descriptor_behler,args_str,error)
      case(DT_DISTANCE_2B)
         call initialise(this%descriptor_distance_2b,args_str,error)
      case(DT_COORDINATION)
         call initialise(this%descriptor_coordination,args_str,error)
      case(DT_ANGLE_3B)
         call initialise(this%descriptor_angle_3b,args_str,error)
      case(DT_CO_ANGLE_3B)
         call initialise(this%descriptor_co_angle_3b,args_str,error)
      case(DT_CO_DISTANCE_2B)
         call initialise(this%descriptor_co_distance_2b,args_str,error)
      case(DT_COSNX)
         call initialise(this%descriptor_cosnx,args_str,error)
      case(DT_TRIHIS)
         call initialise(this%descriptor_trihis,args_str,error)
      case(DT_WATER_MONOMER)
         call initialise(this%descriptor_water_monomer,args_str,error)
      case(DT_WATER_DIMER)
         call initialise(this%descriptor_water_dimer,args_str,error)
      case(DT_A2_DIMER)
         call initialise(this%descriptor_A2_dimer,args_str,error)
      case(DT_AB_DIMER)
         call initialise(this%descriptor_AB_dimer,args_str,error)
      case(DT_BOND_REAL_SPACE)
         call initialise(this%descriptor_bond_real_space,args_str,error)
      case(DT_ATOM_REAL_SPACE)
         call initialise(this%descriptor_atom_real_space,args_str,error)
      case(DT_POWER_SO3)
         call initialise(this%descriptor_power_so3,args_str,error)
      case(DT_POWER_SO4)
         call initialise(this%descriptor_power_so4,args_str,error)
      case(DT_SOAP)
         call initialise(this%descriptor_soap,args_str,error)
      case(DT_AN_MONOMER)
         call initialise(this%descriptor_AN_monomer,args_str,error)
      case(DT_GENERAL_MONOMER)
         call initialise(this%descriptor_general_monomer,args_str,error)
      case(DT_GENERAL_DIMER)
         call initialise(this%descriptor_general_dimer,args_str,error)
      case(DT_GENERAL_TRIMER)
         call initialise(this%descriptor_general_trimer,args_str,error)
      case(DT_RDF)
         call initialise(this%descriptor_rdf,args_str,error)
      case(DT_AS_DISTANCE_2B)
         call initialise(this%descriptor_as_distance_2b,args_str,error)
      case(DT_MOLECULE_LO_D)
         call initialise(this%descriptor_molecule_lo_d,args_str,error)
      case(DT_ALEX)
         call initialise(this%descriptor_alex,args_str,error)
      case(DT_COM_DIMER)
         call initialise(this%descriptor_com_dimer,args_str,error)
      case(DT_DISTANCE_NB)
         call initialise(this%descriptor_distance_Nb,args_str,error)
      endselect

   endsubroutine descriptor_initialise

   subroutine descriptor_finalise(this,error)
      type(descriptor), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call finalise(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            call finalise(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            call finalise(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            call finalise(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            call finalise(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            call finalise(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            call finalise(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            call finalise(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            call finalise(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            call finalise(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            call finalise(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            call finalise(this%descriptor_water_dimer,error)
         case(DT_A2_dimer)
            call finalise(this%descriptor_A2_dimer,error)
         case(DT_AB_dimer)
            call finalise(this%descriptor_AB_dimer,error)
         case(DT_BOND_REAL_SPACE)
            call finalise(this%descriptor_bond_real_space,error)
         case(DT_ATOM_REAL_SPACE)
            call finalise(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            call finalise(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            call finalise(this%descriptor_power_so4,error)
         case(DT_SOAP)
            call finalise(this%descriptor_soap,error)
         case(DT_GENERAL_MONOMER)
            call finalise(this%descriptor_general_monomer,error)
         case(DT_GENERAL_DIMER)
            call finalise(this%descriptor_general_dimer,error)
         case(DT_GENERAL_TRIMER)
            call finalise(this%descriptor_general_trimer,error)
         case(DT_RDF)
            call finalise(this%descriptor_rdf,error)
         case(DT_AS_DISTANCE_2b)
            call finalise(this%descriptor_as_distance_2b,error)
         case(DT_MOLECULE_LO_D)
            call finalise(this%descriptor_molecule_lo_d,error)
         case(DT_ALEX)
            call finalise(this%descriptor_alex,error)
         case(DT_COM_DIMER)
            call finalise(this%descriptor_com_dimer,error)
         case(DT_DISTANCE_Nb)
            call finalise(this%descriptor_distance_Nb,error)
      endselect

      this%descriptor_type = DT_NONE

   endsubroutine descriptor_finalise
   
   subroutine descriptor_MPI_setup(this,at,mpi,mpi_mask,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      type(MPI_Context), intent(in) :: mpi
      logical, dimension(:), intent(out) :: mpi_mask
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(mpi%active) then
         select case(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_BISPECTRUM_SO3)
            RAISE_ERROR("descriptor_MPI_setup: bispectrum_so3 not MPI ready.", error)
         case(DT_BEHLER)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_DISTANCE_2B)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_COORDINATION)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_ANGLE_3B)
            RAISE_ERROR("descriptor_MPI_setup: angle_3b not MPI ready.", error)
         case(DT_CO_ANGLE_3B)
            RAISE_ERROR("descriptor_MPI_setup: co_angle_3b not MPI ready.", error)
         case(DT_CO_DISTANCE_2B)
            RAISE_ERROR("descriptor_MPI_setup: co_distance_2b not MPI ready.", error)
         case(DT_COSNX)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_TRIHIS)
            RAISE_ERROR("descriptor_MPI_setup: trihis not MPI ready.", error)
         case(DT_WATER_MONOMER)
            call descriptor_water_monomer_dimer_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_WATER_DIMER)
            call descriptor_water_monomer_dimer_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_A2_DIMER)
            RAISE_ERROR("descriptor_MPI_setup: A2_dimer not MPI ready.", error)
         case(DT_AB_DIMER)
            RAISE_ERROR("descriptor_MPI_setup: AB_dimer not MPI ready.", error)
         case(DT_BOND_REAL_SPACE)
            RAISE_ERROR("descriptor_MPI_setup: bond_real_space not MPI ready.", error)
         case(DT_ATOM_REAL_SPACE)
            RAISE_ERROR("descriptor_MPI_setup: atom_real_space not MPI ready.", error)
         case(DT_POWER_SO3)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_POWER_SO4)
            RAISE_ERROR("descriptor_MPI_setup: power_SO4 not MPI ready.", error)
         case(DT_SOAP)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_AN_MONOMER)
            RAISE_ERROR("descriptor_MPI_setup: AN_monomer not MPI ready.", error)
         case(DT_GENERAL_MONOMER)
            call descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
         case(DT_GENERAL_DIMER)
            call descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
         case(DT_GENERAL_TRIMER)
            call descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
         case(DT_RDF)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_AS_DISTANCE_2B)
            RAISE_ERROR("descriptor_MPI_setup: as_distance_2b not MPI ready.", error)
         case(DT_MOLECULE_LO_D)
            RAISE_ERROR("descriptor_MPI_setup: molecule_lo_d not MPI ready.", error)
         case(DT_ALEX)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_COM_DIMER)
            call descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
         case(DT_DISTANCE_NB)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case default
            RAISE_ERROR("descriptor_MPI_setup: descriptor type "//this%descriptor_type//" not recognised.",error)
         endselect
      else
         mpi_mask = .true.
      endif

   endsubroutine descriptor_MPI_setup

   subroutine descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
      type(atoms), intent(in) :: at
      type(MPI_Context), intent(in) :: mpi
      logical, dimension(:), intent(out) :: mpi_mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)
      
      mpi_mask = .false.
      do i = 1, at%N
         if( mod(i-1, mpi%n_procs) == mpi%my_proc ) mpi_mask(i) = .true.
      enddo

   endsubroutine descriptor_atomic_MPI_setup

   subroutine descriptor_water_monomer_dimer_MPI_setup(at,mpi,mpi_mask,error)
      type(atoms), intent(in) :: at
      type(MPI_Context), intent(in) :: mpi
      logical, dimension(:), intent(out) :: mpi_mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)
      
      mpi_mask = .false.
      do i = 1, at%N
         if( at%Z(i) == 8 .and. mod(i-1, mpi%n_procs) == mpi%my_proc ) mpi_mask(i) = .true.
      enddo

   endsubroutine descriptor_water_monomer_dimer_MPI_setup

   subroutine descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      type(MPI_Context), intent(in) :: mpi
      logical, dimension(:), intent(out) :: mpi_mask
      integer, optional, intent(out) :: error

      integer, dimension(:,:), allocatable :: monomer_index
      logical, dimension(at%N) :: associated_to_monomer
      integer :: n_monomer

      integer :: i

      INIT_ERROR(error)

      associated_to_monomer = .false.
      select case(this%descriptor_type)
      case(DT_GENERAL_MONOMER)
         call find_general_monomer(at,monomer_index,this%descriptor_general_monomer%signature,associated_to_monomer,this%descriptor_general_monomer%cutoff,this%descriptor_general_monomer%atom_ordercheck,error)
      case(DT_GENERAL_DIMER)
         call find_general_monomer(at,monomer_index,this%descriptor_general_dimer%signature_one,associated_to_monomer,this%descriptor_general_dimer%monomer_one_cutoff,this%descriptor_general_dimer%atom_ordercheck,error)
      case(DT_GENERAL_TRIMER)
         call find_general_monomer(at,monomer_index,this%descriptor_general_trimer%signature_one,associated_to_monomer,this%descriptor_general_trimer%monomer_one_cutoff,this%descriptor_general_trimer%atom_ordercheck,error)
      case(DT_COM_DIMER)
         call find_general_monomer(at,monomer_index,this%descriptor_com_dimer%signature_one,associated_to_monomer,this%descriptor_com_dimer%monomer_one_cutoff,this%descriptor_com_dimer%atom_ordercheck,error)
      case default
         RAISE_ERROR("descriptor_general_monomer_nmer_MPI_setup: descriptor type "//this%descriptor_type//" not recognised.",error)
      endselect

      n_monomer = size(monomer_index,2)

      mpi_mask = .false.
      do i = 1, n_monomer ! for dimer, trimer this is the first monomer (signature_one)
         if( mod(i-1, mpi%n_procs) == mpi%my_proc ) then
            mpi_mask(monomer_index(:,i)) = .true.
         endif
      enddo

      deallocate(monomer_index)

   endsubroutine descriptor_general_monomer_nmer_MPI_setup

   subroutine descriptor_data_finalise(this,error)
      type(descriptor_data), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(allocated(this%x)) then
         do i = 1, size(this%x)
            if(allocated(this%x(i)%data)) deallocate(this%x(i)%data)
            if(allocated(this%x(i)%grad_data)) deallocate(this%x(i)%grad_data)
            if(allocated(this%x(i)%ci)) deallocate(this%x(i)%ci)
            if(allocated(this%x(i)%ii)) deallocate(this%x(i)%ii)
            if(allocated(this%x(i)%pos)) deallocate(this%x(i)%pos)
            if(allocated(this%x(i)%has_grad_data)) deallocate(this%x(i)%has_grad_data)
            if(allocated(this%x(i)%grad_covariance_cutoff)) deallocate(this%x(i)%grad_covariance_cutoff)
         enddo
         deallocate(this%x)
      endif

   endsubroutine descriptor_data_finalise

   subroutine RadialFunction_initialise(this,n_max,cutoff, min_cutoff,error)
      type(RadialFunction_type), intent(inout) :: this
      integer, intent(in) :: n_max
      real(dp), intent(in) :: cutoff, min_cutoff
      integer, optional, intent(out) :: error

      real(dp), dimension(:,:), allocatable :: S, vS
      real(dp), dimension(:), allocatable :: eS
      integer :: i, j

      INIT_ERROR(error)

      call finalise(this)

      this%n_max = n_max
      this%cutoff = cutoff
      this%min_cutoff = min_cutoff

      allocate(this%RadialTransform(this%n_max,this%n_max),this%NormFunction(this%n_max))
      allocate(S(this%n_max,this%n_max), vS(this%n_max,this%n_max), eS(this%n_max))

      do i = 1, this%n_max
         this%NormFunction(i) = sqrt(this%cutoff**(2.0_dp*i+5.0_dp)/(2.0_dp*i+5.0_dp))
         do j = 1, this%n_max
            S(j,i) = sqrt((2.0_dp*i+5)*(2.0_dp*j+5))/(i+j+5.0_dp)
         enddo
      enddo

      call diagonalise(S,eS,vS)
      this%RadialTransform = matmul(matmul(vS,diag(1.0_dp/sqrt(eS))),transpose(vS))

      if(allocated(S)) deallocate(S)
      if(allocated(vS)) deallocate(vS)
      if(allocated(eS)) deallocate(eS)

      this%initialised = .true.

   endsubroutine RadialFunction_initialise

   subroutine RadialFunction_finalise(this,error)
      type(RadialFunction_type), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%n_max = 0

      if(allocated(this%RadialTransform)) deallocate(this%RadialTransform)
      if(allocated(this%NormFunction)) deallocate(this%NormFunction)

      this%initialised = .false.

   endsubroutine RadialFunction_finalise

   subroutine cplx_2d_array1_finalise(this)
      type(cplx_2d), dimension(:), allocatable, intent(inout) :: this
      integer :: j

      if(allocated(this)) then
         do j = lbound(this,1), ubound(this,1)
            if(allocated(this(j)%mm)) deallocate(this(j)%mm)
         enddo
         deallocate(this)
      endif
   endsubroutine cplx_2d_array1_finalise

   subroutine cplx_3d_array2_finalise(this)
      type(cplx_3d), dimension(:,:), allocatable, intent(inout) :: this
      integer :: i, j

      if(allocated(this)) then
         do j = lbound(this,2), ubound(this,2)
            do i = lbound(this,1), ubound(this,1)
               if(allocated(this(i,j)%mm)) deallocate(this(i,j)%mm)
            enddo
         enddo
         deallocate(this)
      endif

   endsubroutine cplx_3d_array2_finalise

   subroutine fourier_SO4_calc(this,at,i,U,dU,args_str,error)
      type(fourier_SO4_type), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(in) :: i
      type(cplx_2d), dimension(:), allocatable, intent(inout) :: U
      type(cplx_3d), dimension(:,:), allocatable, intent(inout), optional :: dU
      integer, optional, intent(out) :: error
      character(len=*), intent(in), optional :: args_str 

      complex(dp), dimension(:,:), allocatable :: Uc, Up
      complex(dp), dimension(:,:,:), allocatable :: dUc, dUp
      complex(dp) :: z0_pls_Iz, z0_min_Iz, x_pls_Iy, x_min_Iy
      complex(dp), dimension(3) :: dz0_pls_Iz, dz0_min_Iz, dx_pls_Iy, dx_min_Iy
      real(dp), dimension(3) :: diff, u_ij, dfcut, dz0, dr0
      real(dp) :: r0, r, fcut, z0, theta0
      integer :: n, n_i, ji, j, m1, m2
      integer, dimension(116) :: species_map

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR('fourier_SO4_calc: object not initialised',error)
      endif

      species_map = 0
      do j = 1, size(this%species_Z)
         if(this%species_Z(j) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(j)) = j
         endif
      enddo


      if(allocated(U)) then
         if(lbound(U,1) /= 0 .or. ubound(U,1) /= this%j_max) call finalise(U)
      endif

      if(.not.allocated(U)) then
         allocate( U(0:this%j_max) )
         do j = 0, this%j_max
            allocate( U(j)%mm(-j:j,-j:j) )
            U(j)%mm = CPLX_ZERO
         enddo
      endif

      do j = 0, this%j_max
         U(j)%mm = CPLX_ZERO
         do m1 = -j, j, 2
            U(j)%mm(m1,m1) = CPLX_ONE
         enddo
      enddo

      allocate( Uc(-this%j_max:this%j_max, -this%j_max:this%j_max), &
        Up(-this%j_max:this%j_max, -this%j_max:this%j_max) )

      Uc = CPLX_ZERO
      Up = CPLX_ZERO

      if(present(dU)) then
         if(allocated(dU)) call finalise(dU)

         ! dU is not allocated, allocate and zero it
         allocate( dU(0:this%j_max,0:n_neighbours(at,i,max_dist=this%cutoff)) )
         do j = 0, this%j_max
            allocate( dU(j,0)%mm(3,-j:j,-j:j) )
            dU(j,0)%mm = CPLX_ZERO
         enddo

         allocate( dUc(3,-this%j_max:this%j_max, -this%j_max:this%j_max), &
            dUp(3,-this%j_max:this%j_max, -this%j_max:this%j_max) )
         dUc = CPLX_ZERO
         dUp = CPLX_ZERO
      endif

      n_i = 0
      do n = 1, n_neighbours(at,i)
         ji = neighbour(at, i, n, distance=r, diff=diff, cosines=u_ij)
         if( r >= this%cutoff ) cycle

         n_i = n_i + 1

         theta0 = r / this%z0
         z0 = r / tan( theta0 )
         r0 = sin( theta0 ) / r

         z0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) ) * r0
         z0_min_Iz = ( z0 - CPLX_IMAG*diff(3) ) * r0
         x_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) ) * r0
         x_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) ) * r0

         fcut = cos_cutoff_function(r,this%cutoff) * this%w(species_map(at%Z(ji)))

         U(0)%mm(0,0) = U(0)%mm(0,0) + fcut
         Up(0:0,0:0) = CPLX_ONE

         if(present(dU)) then

            dfcut = -dcos_cutoff_function(r,this%cutoff)*u_ij * this%w(species_map(at%Z(ji)))
            dz0 = ( 1.0_dp / tan( theta0 ) - theta0 / sin(theta0)**2 ) * u_ij
            dr0 = ( cos( theta0 ) / (r*this%z0) - r0 / r ) * u_ij

            dz0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) )*dr0 + dz0*r0
            dz0_pls_Iz(3) = dz0_pls_Iz(3) + CPLX_IMAG*r0
               
            dz0_min_Iz = ( z0 - CPLX_IMAG*diff(3) )*dr0 + dz0*r0
            dz0_min_Iz(3) = dz0_min_Iz(3) - CPLX_IMAG*r0
               
            dx_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) )*dr0
            dx_pls_Iy(1) = dx_pls_Iy(1) + r0
            dx_pls_Iy(2) = dx_pls_Iy(2) + CPLX_IMAG*r0
               
            dx_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) )*dr0
            dx_min_Iy(1) = dx_min_Iy(1) + r0
            dx_min_Iy(2) = dx_min_Iy(2) - CPLX_IMAG*r0

            dUc = CPLX_ZERO
            dUp = CPLX_ZERO

            dU(0,0)%mm(:,0,0) = dU(0,0)%mm(:,0,0) + dfcut*CPLX_ONE

            allocate( dU(0,n_i)%mm(3,-0:0,-0:0) )

            dU(0,n_i)%mm(:,0,0) = - dfcut*CPLX_ONE
         endif

         do j = 1, this%j_max
            Uc(-j:j,-j:j) = CPLX_ZERO
            if(present(dU)) then
               dUc(:,-j:j,-j:j) = CPLX_ZERO
               allocate( dU(j,n_i)%mm(3,-j:j,-j:j) )
               dU(j,n_i)%mm = CPLX_ZERO
            endif

            do m1 = -j, j-2, 2
               do m2 = -j, j, 2
                  if( (j-m2) /= 0 ) then
                     Uc(m2,m1) = Uc(m2,m1) + &
                     sqrt( real(j-m2,dp)/real(j-m1,dp) ) * z0_pls_Iz * Up(m2+1,m1+1)

                     if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) + &
                        sqrt( real(j-m2,dp)/real(j-m1,dp) ) * &
                        ( dz0_pls_Iz * Up(m2+1,m1+1) + z0_pls_Iz * dUp(:,m2+1,m1+1) )
                  endif

                  if( (j+m2) /= 0 ) then
                     Uc(m2,m1) = Uc(m2,m1) - &
                     CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * x_min_Iy * Up(m2-1,m1+1)

                     if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) - &
                        CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * &
                        ( dx_min_Iy * Up(m2-1,m1+1) + x_min_Iy * dUp(:,m2-1,m1+1) )

                  endif
               enddo
            enddo

            m1 = j
            do m2 = -j, j, 2
               if( (j+m2) /= 0 ) then
                  Uc(m2,m1) = Uc(m2,m1) + &
                  sqrt( real(j+m2,dp)/real(j+m1,dp) ) * z0_min_Iz * Up(m2-1,m1-1)

                  if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) + &
                     sqrt( real(j+m2,dp)/real(j+m1,dp) ) * &
                     ( dz0_min_Iz * Up(m2-1,m1-1) + z0_min_Iz * dUp(:,m2-1,m1-1) )
               endif

               if( (j-m2) /= 0 ) then
                  Uc(m2,m1) = Uc(m2,m1) - &
                  CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * x_pls_Iy * Up(m2+1,m1-1)

                  if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) - &
                     CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * &
                     ( dx_pls_Iy * Up(m2+1,m1-1) + x_pls_Iy * dUp(:,m2+1,m1-1) )
               endif
            enddo

            U(j)%mm = U(j)%mm + Uc(-j:j,-j:j) * fcut
            Up(-j:j,-j:j) = Uc(-j:j,-j:j)
            if(present(dU)) then
               dUp(:,-j:j,-j:j) = dUc(:,-j:j,-j:j)
               dU(j,0)%mm = dU(j,0)%mm - dUc(:,-j:j,-j:j) * fcut
               dU(j,n_i)%mm = dU(j,n_i)%mm + dUc(:,-j:j,-j:j) * fcut
               do m1 = -j, j, 2
                  do m2 = -j, j, 2
                     dU(j,0)%mm(:,m2,m1) = dU(j,0)%mm(:,m2,m1) &
                        + Uc(m2,m1) * dfcut
                     dU(j,n_i)%mm(:,m2,m1) = dU(j,n_i)%mm(:,m2,m1) &
                        - Uc(m2,m1) * dfcut
                  enddo
               enddo
            endif

         enddo ! j
      enddo ! n

      if(allocated(Up)) deallocate(Up)
      if(allocated(Uc)) deallocate(Uc)
      if(allocated(dUp)) deallocate(dUp)
      if(allocated(dUc)) deallocate(dUc)

   endsubroutine fourier_SO4_calc

   subroutine fourier_so4_initialise(this,args_str,error)
      type(fourier_SO4_type), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '2.75', this%cutoff, help_string="Cutoff for SO4 bispectrum")
      call param_register(params, 'z0_ratio', '0.0', this%z0_ratio, help_string="Ratio of radius of 4D projection sphere times PI and the cutoff.")
      call param_register(params, 'j_max', '4', this%j_max, help_string="Max of expansion of bispectrum, i.e. resulution")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='fourier_so4_initialise args_str')) then
         RAISE_ERROR("fourier_so4_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='fourier_so4_initialise args_str')) then
         RAISE_ERROR("fourier_so4_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%z0 = max(1.0_dp,this%z0_ratio) * this%cutoff/(PI-0.02_dp)
      
      this%initialised = .true.


   endsubroutine fourier_so4_initialise

   subroutine fourier_so4_finalise(this,error)
      type(fourier_so4_type), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%j_max = 0
      this%z0_ratio = 0.0_dp
      this%z0 = 0.0_dp
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      this%initialised = .false.

   endsubroutine fourier_so4_finalise

   subroutine bispectrum_so4_initialise(this,args_str,error)
      type(bispectrum_so4), intent(inout), target :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      call initialise(this%fourier_SO4,args_str,error)

      this%cutoff => this%fourier_SO4%cutoff
      this%z0_ratio => this%fourier_SO4%z0_ratio
      this%z0 => this%fourier_SO4%z0
      this%j_max => this%fourier_SO4%j_max
      this%Z => this%fourier_SO4%Z
      this%cutoff => this%fourier_SO4%cutoff
      this%species_Z => this%fourier_SO4%species_Z
      this%w => this%fourier_SO4%w
      
      this%initialised = .true.

   endsubroutine bispectrum_so4_initialise

   subroutine bispectrum_so4_finalise(this,error)
      type(bispectrum_so4), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      call finalise(this%fourier_SO4,error)

      this%cutoff => null()
      this%z0_ratio => null()
      this%z0 => null()
      this%j_max => null()
      this%Z => null()
      this%cutoff => null()
      this%species_Z => null()
      this%w => null()

      this%initialised = .false.

   endsubroutine bispectrum_so4_finalise

   subroutine bispectrum_so3_initialise(this,args_str,error)
      type(bispectrum_so3), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)

      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for bispectrum_so3-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in bispectrum_so3-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for bispectrum_so3-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for bispectrum_so3-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bispectrum_so3_initialise args_str')) then
         RAISE_ERROR("bispectrum_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bispectrum_so3_initialise args_str')) then
         RAISE_ERROR("bispectrum_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

      call print('Dimensions: '//bispectrum_so3_dimensions(this,error))

   endsubroutine bispectrum_so3_initialise

   subroutine bispectrum_so3_finalise(this,error)
      type(bispectrum_so3), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine bispectrum_so3_finalise

   subroutine behler_initialise(this,args_str,error)
      type(behler), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      !call param_register(params, 'behler_cutoff', '2.75', this%cutoff, help_string="Cutoff for Behler-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='behler_initialise args_str')) then
         RAISE_ERROR("behler_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif

      call finalise(params)

      this%n_g2 = 8
      this%n_g3 = 43

      allocate(this%g2(this%n_g2), this%g3(this%n_g3))

      this%g2(1)%eta = 0.001_dp / BOHR**2; this%g2(1)%rs = 0.000_dp * BOHR; this%g2(1)%rc = 11.338_dp * BOHR
      this%g2(2)%eta = 0.010_dp / BOHR**2; this%g2(2)%rs = 0.000_dp * BOHR; this%g2(2)%rc = 11.338_dp * BOHR
      this%g2(3)%eta = 0.020_dp / BOHR**2; this%g2(3)%rs = 0.000_dp * BOHR; this%g2(3)%rc = 11.338_dp * BOHR
      this%g2(4)%eta = 0.035_dp / BOHR**2; this%g2(4)%rs = 0.000_dp * BOHR; this%g2(4)%rc = 11.338_dp * BOHR
      this%g2(5)%eta = 0.060_dp / BOHR**2; this%g2(5)%rs = 0.000_dp * BOHR; this%g2(5)%rc = 11.338_dp * BOHR
      this%g2(6)%eta = 0.100_dp / BOHR**2; this%g2(6)%rs = 0.000_dp * BOHR; this%g2(6)%rc = 11.338_dp * BOHR
      this%g2(7)%eta = 0.200_dp / BOHR**2; this%g2(7)%rs = 0.000_dp * BOHR; this%g2(7)%rc = 11.338_dp * BOHR
      this%g2(8)%eta = 0.400_dp / BOHR**2; this%g2(8)%rs = 0.000_dp * BOHR; this%g2(8)%rc = 11.338_dp * BOHR

      this%g3( 1)%eta = 0.0001_dp / BOHR**2; this%g3( 1)%lambda = -1.000_dp; this%g3( 1)%zeta =  1.000_dp; this%g3( 1)%rc = 11.338_dp * BOHR
      this%g3( 2)%eta = 0.0001_dp / BOHR**2; this%g3( 2)%lambda =  1.000_dp; this%g3( 2)%zeta =  1.000_dp; this%g3( 2)%rc = 11.338_dp * BOHR
      this%g3( 3)%eta = 0.0001_dp / BOHR**2; this%g3( 3)%lambda = -1.000_dp; this%g3( 3)%zeta =  2.000_dp; this%g3( 3)%rc = 11.338_dp * BOHR
      this%g3( 4)%eta = 0.0001_dp / BOHR**2; this%g3( 4)%lambda =  1.000_dp; this%g3( 4)%zeta =  2.000_dp; this%g3( 4)%rc = 11.338_dp * BOHR
      this%g3( 5)%eta = 0.0030_dp / BOHR**2; this%g3( 5)%lambda = -1.000_dp; this%g3( 5)%zeta =  1.000_dp; this%g3( 5)%rc = 11.338_dp * BOHR
      this%g3( 6)%eta = 0.0030_dp / BOHR**2; this%g3( 6)%lambda =  1.000_dp; this%g3( 6)%zeta =  1.000_dp; this%g3( 6)%rc = 11.338_dp * BOHR
      this%g3( 7)%eta = 0.0030_dp / BOHR**2; this%g3( 7)%lambda = -1.000_dp; this%g3( 7)%zeta =  2.000_dp; this%g3( 7)%rc = 11.338_dp * BOHR
      this%g3( 8)%eta = 0.0030_dp / BOHR**2; this%g3( 8)%lambda =  1.000_dp; this%g3( 8)%zeta =  2.000_dp; this%g3( 8)%rc = 11.338_dp * BOHR
      this%g3( 9)%eta = 0.0080_dp / BOHR**2; this%g3( 9)%lambda = -1.000_dp; this%g3( 9)%zeta =  1.000_dp; this%g3( 9)%rc = 11.338_dp * BOHR
      this%g3(10)%eta = 0.0080_dp / BOHR**2; this%g3(10)%lambda =  1.000_dp; this%g3(10)%zeta =  1.000_dp; this%g3(10)%rc = 11.338_dp * BOHR
      this%g3(11)%eta = 0.0080_dp / BOHR**2; this%g3(11)%lambda = -1.000_dp; this%g3(11)%zeta =  2.000_dp; this%g3(11)%rc = 11.338_dp * BOHR
      this%g3(12)%eta = 0.0080_dp / BOHR**2; this%g3(12)%lambda =  1.000_dp; this%g3(12)%zeta =  2.000_dp; this%g3(12)%rc = 11.338_dp * BOHR
      this%g3(13)%eta = 0.0150_dp / BOHR**2; this%g3(13)%lambda = -1.000_dp; this%g3(13)%zeta =  1.000_dp; this%g3(13)%rc = 11.338_dp * BOHR
      this%g3(14)%eta = 0.0150_dp / BOHR**2; this%g3(14)%lambda =  1.000_dp; this%g3(14)%zeta =  1.000_dp; this%g3(14)%rc = 11.338_dp * BOHR
      this%g3(15)%eta = 0.0150_dp / BOHR**2; this%g3(15)%lambda = -1.000_dp; this%g3(15)%zeta =  2.000_dp; this%g3(15)%rc = 11.338_dp * BOHR
      this%g3(16)%eta = 0.0150_dp / BOHR**2; this%g3(16)%lambda =  1.000_dp; this%g3(16)%zeta =  2.000_dp; this%g3(16)%rc = 11.338_dp * BOHR
      this%g3(17)%eta = 0.0150_dp / BOHR**2; this%g3(17)%lambda = -1.000_dp; this%g3(17)%zeta =  4.000_dp; this%g3(17)%rc = 11.338_dp * BOHR
      this%g3(18)%eta = 0.0150_dp / BOHR**2; this%g3(18)%lambda =  1.000_dp; this%g3(18)%zeta =  4.000_dp; this%g3(18)%rc = 11.338_dp * BOHR
      this%g3(19)%eta = 0.0150_dp / BOHR**2; this%g3(19)%lambda = -1.000_dp; this%g3(19)%zeta = 16.000_dp; this%g3(19)%rc = 11.338_dp * BOHR
      this%g3(20)%eta = 0.0150_dp / BOHR**2; this%g3(20)%lambda =  1.000_dp; this%g3(20)%zeta = 16.000_dp; this%g3(20)%rc = 11.338_dp * BOHR
      this%g3(21)%eta = 0.0250_dp / BOHR**2; this%g3(21)%lambda = -1.000_dp; this%g3(21)%zeta =  1.000_dp; this%g3(21)%rc = 11.338_dp * BOHR
      this%g3(22)%eta = 0.0250_dp / BOHR**2; this%g3(22)%lambda =  1.000_dp; this%g3(22)%zeta =  1.000_dp; this%g3(22)%rc = 11.338_dp * BOHR
      this%g3(23)%eta = 0.0250_dp / BOHR**2; this%g3(23)%lambda = -1.000_dp; this%g3(23)%zeta =  2.000_dp; this%g3(23)%rc = 11.338_dp * BOHR
      this%g3(24)%eta = 0.0250_dp / BOHR**2; this%g3(24)%lambda =  1.000_dp; this%g3(24)%zeta =  2.000_dp; this%g3(24)%rc = 11.338_dp * BOHR
      this%g3(25)%eta = 0.0250_dp / BOHR**2; this%g3(25)%lambda = -1.000_dp; this%g3(25)%zeta =  4.000_dp; this%g3(25)%rc = 11.338_dp * BOHR
      this%g3(26)%eta = 0.0250_dp / BOHR**2; this%g3(26)%lambda =  1.000_dp; this%g3(26)%zeta =  4.000_dp; this%g3(26)%rc = 11.338_dp * BOHR
      this%g3(27)%eta = 0.0250_dp / BOHR**2; this%g3(27)%lambda = -1.000_dp; this%g3(27)%zeta = 16.000_dp; this%g3(27)%rc = 11.338_dp * BOHR
      this%g3(28)%eta = 0.0250_dp / BOHR**2; this%g3(28)%lambda =  1.000_dp; this%g3(28)%zeta = 16.000_dp; this%g3(28)%rc = 11.338_dp * BOHR
      this%g3(29)%eta = 0.0450_dp / BOHR**2; this%g3(29)%lambda = -1.000_dp; this%g3(29)%zeta =  1.000_dp; this%g3(29)%rc = 11.338_dp * BOHR
      this%g3(30)%eta = 0.0450_dp / BOHR**2; this%g3(30)%lambda =  1.000_dp; this%g3(30)%zeta =  1.000_dp; this%g3(30)%rc = 11.338_dp * BOHR
      this%g3(31)%eta = 0.0450_dp / BOHR**2; this%g3(31)%lambda = -1.000_dp; this%g3(31)%zeta =  2.000_dp; this%g3(31)%rc = 11.338_dp * BOHR
      this%g3(32)%eta = 0.0450_dp / BOHR**2; this%g3(32)%lambda =  1.000_dp; this%g3(32)%zeta =  2.000_dp; this%g3(32)%rc = 11.338_dp * BOHR
      this%g3(33)%eta = 0.0450_dp / BOHR**2; this%g3(33)%lambda = -1.000_dp; this%g3(33)%zeta =  4.000_dp; this%g3(33)%rc = 11.338_dp * BOHR
      this%g3(34)%eta = 0.0450_dp / BOHR**2; this%g3(34)%lambda =  1.000_dp; this%g3(34)%zeta =  4.000_dp; this%g3(34)%rc = 11.338_dp * BOHR
      this%g3(35)%eta = 0.0450_dp / BOHR**2; this%g3(35)%lambda = -1.000_dp; this%g3(35)%zeta = 16.000_dp; this%g3(35)%rc = 11.338_dp * BOHR
      this%g3(36)%eta = 0.0450_dp / BOHR**2; this%g3(36)%lambda =  1.000_dp; this%g3(36)%zeta = 16.000_dp; this%g3(36)%rc = 11.338_dp * BOHR
      this%g3(37)%eta = 0.0800_dp / BOHR**2; this%g3(37)%lambda = -1.000_dp; this%g3(37)%zeta =  1.000_dp; this%g3(37)%rc = 11.338_dp * BOHR
      this%g3(38)%eta = 0.0800_dp / BOHR**2; this%g3(38)%lambda =  1.000_dp; this%g3(38)%zeta =  1.000_dp; this%g3(38)%rc = 11.338_dp * BOHR
      this%g3(39)%eta = 0.0800_dp / BOHR**2; this%g3(39)%lambda = -1.000_dp; this%g3(39)%zeta =  2.000_dp; this%g3(39)%rc = 11.338_dp * BOHR
      this%g3(40)%eta = 0.0800_dp / BOHR**2; this%g3(40)%lambda =  1.000_dp; this%g3(40)%zeta =  2.000_dp; this%g3(40)%rc = 11.338_dp * BOHR
      this%g3(41)%eta = 0.0800_dp / BOHR**2; this%g3(41)%lambda = -1.000_dp; this%g3(41)%zeta =  4.000_dp; this%g3(41)%rc = 11.338_dp * BOHR
      this%g3(42)%eta = 0.0800_dp / BOHR**2; this%g3(42)%lambda =  1.000_dp; this%g3(42)%zeta =  4.000_dp; this%g3(42)%rc = 11.338_dp * BOHR
      this%g3(43)%eta = 0.0800_dp / BOHR**2; this%g3(43)%lambda =  1.000_dp; this%g3(43)%zeta = 16.000_dp; this%g3(43)%rc = 11.338_dp * BOHR

      this%cutoff = 11.338_dp * BOHR

      this%initialised = .true.

   endsubroutine behler_initialise

   subroutine behler_finalise(this,error)
      type(behler), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%n_g2 = 0
      this%n_g3 = 0

      if(allocated(this%g2)) deallocate(this%g2)
      if(allocated(this%g3)) deallocate(this%g3)

      this%initialised = .false.

   endsubroutine behler_finalise

   subroutine distance_2b_initialise(this,args_str,error)
      type(distance_2b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      logical :: has_resid_name

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for distance_2b-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.5', this%cutoff_transition_width, help_string="Transition width of cutoff for distance_2b-type descriptors")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atom type #1 in bond")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atom type #2 in bond")
      call param_register(params, 'resid_name', '', this%resid_name, has_value_target=has_resid_name, help_string="Name of an integer property in the atoms object giving the residue id of the molecule to which the atom belongs.")
      call param_register(params, 'only_intra', 'F', this%only_intra, help_string="Only calculate INTRAmolecular pairs with equal residue ids (bonds)")
      call param_register(params, 'only_inter', 'F', this%only_inter, help_string="Only apply to INTERmolecular pairs with different residue ids (non-bonded)")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='distance_2b_initialise args_str')) then
         RAISE_ERROR("distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if (this%only_intra .and. this%only_inter) then
         RAISE_ERROR("distance_2b_initialise: cannot specify both only_inter AND only_intra", error)
      end if
      if ((this%only_intra .or. this%only_inter) .and. (.not. has_resid_name)) then
         RAISE_ERROR("distance_2b_initialise: only_intra and only_inter require resid_name to be given as well", error)
      end if

      this%initialised = .true.

   endsubroutine distance_2b_initialise

   subroutine distance_2b_finalise(this,error)
      type(distance_2b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.5_dp
      this%Z1 = 0
      this%Z2 = 0

      this%resid_name = ''
      this%only_intra = .false.
      this%only_inter = .false.

      this%initialised = .false.

   endsubroutine distance_2b_finalise

   subroutine coordination_initialise(this,args_str,error)
      type(coordination), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)
      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for coordination-type descriptors")
      call param_register(params, 'transition_width', '0.20', this%transition_width, help_string="Width of transition region from 1 to 0")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='coordination_initialise args_str')) then
         RAISE_ERROR("coordination_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine coordination_initialise

   subroutine coordination_finalise(this,error)
      type(coordination), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%transition_width = 0.0_dp
      this%Z = 0

      this%initialised = .false.

   endsubroutine coordination_finalise

   subroutine angle_3b_initialise(this,args_str,error)
      type(angle_3b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for angle_3b-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atomic number of neighbour #1")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atomic number of neighbour #2")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='angle_3b_initialise args_str')) then
         RAISE_ERROR("angle_3b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine angle_3b_initialise

   subroutine angle_3b_finalise(this,error)
      type(angle_3b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%Z = 0
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine angle_3b_finalise

   subroutine co_angle_3b_initialise(this,args_str,error)
      type(co_angle_3b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for co_angle_3b-type descriptors")
      call param_register(params, 'coordination_cutoff', '0.00', this%coordination_cutoff, help_string="Cutoff for coordination function in co_angle_3b-type descriptors")
      call param_register(params, 'coordination_transition_width', '0.00', this%coordination_transition_width, help_string="Transition width for co_angle_3b-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atomic number of neighbour #1")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atomic number of neighbour #2")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='co_angle_3b_initialise args_str')) then
         RAISE_ERROR("co_angle_3b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine co_angle_3b_initialise

   subroutine co_angle_3b_finalise(this,error)
      type(co_angle_3b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%coordination_cutoff = 0.0_dp
      this%coordination_transition_width = 0.0_dp
      this%Z = 0
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine co_angle_3b_finalise

   subroutine co_distance_2b_initialise(this,args_str,error)
      type(co_distance_2b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for co_distance_2b-type descriptors")
      call param_register(params, 'transition_width', '0.50', this%transition_width, help_string="Transition width of cutoff for co_distance_2b-type descriptors")
      call param_register(params, 'coordination_cutoff', '0.00', this%coordination_cutoff, help_string="Cutoff for coordination function in co_distance_2b-type descriptors")
      call param_register(params, 'coordination_transition_width', '0.00', this%coordination_transition_width, help_string="Transition width for co_distance_2b-type descriptors")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atom type #1 in bond")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atom type #2 in bond")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='co_distance_2b_initialise args_str')) then
         RAISE_ERROR("co_distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine co_distance_2b_initialise

   subroutine co_distance_2b_finalise(this,error)
      type(co_distance_2b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%coordination_cutoff = 0.0_dp
      this%coordination_transition_width = 0.0_dp
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine co_distance_2b_finalise

   subroutine cosnx_initialise(this,args_str,error)
      type(cosnx), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for cosnx-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in cosnx-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for cosnx-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for cosnx-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='cosnx_initialise args_str')) then
         RAISE_ERROR("cosnx_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='cosnx_initialise args_str')) then
         RAISE_ERROR("cosnx_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

   endsubroutine cosnx_initialise

   subroutine cosnx_finalise(this,error)
      type(cosnx), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine cosnx_finalise

   subroutine trihis_initialise(this,args_str,error)
      type(trihis), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      real(dp), dimension(:), allocatable :: gauss_centre1D, gauss_width1D

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for trihis-type descriptors")
      call param_register(params, 'n_gauss', '0', this%n_gauss, help_string="Number of Gaussians for trihis-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='trihis_initialise args_str')) then
         RAISE_ERROR("trihis_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(gauss_centre1D(3*this%n_gauss),gauss_width1D(3*this%n_gauss))
      allocate(this%gauss_centre(3,this%n_gauss),this%gauss_width(3,this%n_gauss))

      call initialise(params)
      call param_register(params, 'trihis_gauss_centre', PARAM_MANDATORY, gauss_centre1D, help_string="Number of Gaussians for trihis-type descriptors")
      call param_register(params, 'trihis_gauss_width', PARAM_MANDATORY, gauss_width1D, help_string="Number of Gaussians for trihis-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='trihis_initialise args_str')) then
         RAISE_ERROR("trihis_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%gauss_centre = reshape(gauss_centre1D,(/3,this%n_gauss/))
      this%gauss_width = reshape(gauss_width1D,(/3,this%n_gauss/))

      deallocate(gauss_centre1D,gauss_width1D)

      this%initialised = .true.

   endsubroutine trihis_initialise

   subroutine trihis_finalise(this,error)
      type(trihis), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%n_gauss = 0

      if(allocated(this%gauss_centre)) deallocate(this%gauss_centre)
      if(allocated(this%gauss_width)) deallocate(this%gauss_width)

      this%initialised = .false.

   endsubroutine trihis_finalise

   subroutine water_monomer_initialise(this,args_str,error)
      type(water_monomer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for water_monomer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='water_monomer_initialise args_str')) then
         RAISE_ERROR("water_monomer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine water_monomer_initialise

   subroutine water_monomer_finalise(this,error)
      type(water_monomer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp

      this%initialised = .false.

   endsubroutine water_monomer_finalise

   subroutine water_dimer_initialise(this,args_str,error)
      type(water_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for water_dimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Width of smooth cutoff region for water_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for water_dimer-type descriptors")
      call param_register(params, 'OHH_ordercheck', 'T', this%OHH_ordercheck, help_string="T: find water molecules. F: use default order OHH")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='water_dimer_initialise args_str')) then
         RAISE_ERROR("water_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine water_dimer_initialise

   subroutine water_dimer_finalise(this,error)
      type(water_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%OHH_ordercheck = .true.

      this%initialised = .false.

   endsubroutine water_dimer_finalise

   subroutine A2_dimer_initialise(this,args_str,error)
      type(A2_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for A2_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for A2_dimer-type descriptors")
      call param_register(params, 'atomic_number', '1', this%atomic_number, help_string="Atomic number in A2_dimer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='A2_dimer_initialise args_str')) then
         RAISE_ERROR("A2_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine A2_dimer_initialise

   subroutine A2_dimer_finalise(this,error)
      type(A2_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%atomic_number = 0

      this%initialised = .false.

   endsubroutine A2_dimer_finalise

   subroutine AB_dimer_initialise(this,args_str,error)
      type(AB_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for AB_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for AB_dimer-type descriptors")
      call param_register(params, 'atomic_number1', '1', this%atomic_number1, help_string="Atomic number of atom 1 in AB_dimer-type descriptors")
      call param_register(params, 'atomic_number2', '9', this%atomic_number2, help_string="Atomic number of atom 2 in AB_dimer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='AB_dimer_initialise args_str')) then
         RAISE_ERROR("AB_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if( this%atomic_number1 == this%atomic_number2 ) then
         RAISE_ERROR("AB_dimer_initialise: AB_dimer_atomic_number1 = AB_dimer_atomic_number2 = "//this%atomic_number1//" which would require addtional permutational symmetries. Use A2_dimer descriptor instead.",error)
      endif

      this%initialised = .true.

   endsubroutine AB_dimer_initialise

   subroutine AB_dimer_finalise(this,error)
      type(AB_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%atomic_number1 = 0
      this%atomic_number2 = 0

      this%initialised = .false.

   endsubroutine AB_dimer_finalise

   subroutine bond_real_space_initialise(this,args_str,error)
      type(bond_real_space), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'bond_cutoff', '0.00', this%bond_cutoff, help_string="Bond cutoff for bond_real_space-type descriptors")
      call param_register(params, 'bond_transition_width', '0.00', this%bond_transition_width, help_string="Bond transition width for bond_real_space-type descriptors")
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Space cutoff for bond_real_space-type descriptors")
      call param_register(params, 'transition_width', '0.00', this%transition_width, help_string="Space transition width for bond_real_space-type descriptors")
      call param_register(params, 'atom_sigma', '0.00', this%atom_sigma, help_string="Atom sigma for bond_real_space-type descriptors")
      call param_register(params, 'max_neighbours', '0', this%max_neighbours, help_string="Maximum number of neighbours")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bond_real_space_initialise args_str')) then
         RAISE_ERROR("bond_real_space_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine bond_real_space_initialise

   subroutine bond_real_space_finalise(this,error)
      type(bond_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%bond_cutoff = 0.0_dp
      this%bond_transition_width = 0.0_dp
      this%cutoff = 0.0_dp
      this%transition_width = 0.0_dp
      this%atom_sigma = 0.0_dp
      this%max_neighbours = 0

      this%initialised = .false.

   endsubroutine bond_real_space_finalise

   subroutine atom_real_space_initialise(this,args_str,error)
      type(atom_real_space), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Space cutoff for atom_real_space-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.00', this%cutoff_transition_width, help_string="Space transition width for atom_real_space-type descriptors")
      call param_register(params, 'l_max', '0', this%l_max, help_string="Cutoff for spherical harmonics expansion")
      call param_register(params, 'alpha', '1.0', this%alpha, help_string="Width of atomic Gaussians")
      call param_register(params, 'zeta', '1.0', this%zeta, help_string="Exponent of covariance function")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='atom_real_space_initialise args_str')) then
         RAISE_ERROR("atom_real_space_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine atom_real_space_initialise

   subroutine atom_real_space_finalise(this,error)
      type(atom_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%l_max = 0
      this%alpha = 0.0_dp
      this%zeta = 0.0_dp

      this%initialised = .false.

   endsubroutine atom_real_space_finalise

   subroutine power_so3_initialise(this,args_str,error)
      type(power_so3), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for power_so3-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in power_so3-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for power_so3-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for power_so3-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='power_so3_initialise args_str')) then
         RAISE_ERROR("power_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='power_so3_initialise args_str')) then
         RAISE_ERROR("power_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

   endsubroutine power_so3_initialise

   subroutine power_so3_finalise(this,error)
      type(power_so3), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine power_so3_finalise

   subroutine power_so4_initialise(this,args_str,error)
      type(power_so4), intent(inout), target :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      call initialise(this%fourier_SO4,args_str,error)

      this%cutoff => this%fourier_SO4%cutoff
      this%z0_ratio => this%fourier_SO4%z0_ratio
      this%z0 => this%fourier_SO4%z0
      this%j_max => this%fourier_SO4%j_max
      this%Z => this%fourier_SO4%Z
      this%cutoff => this%fourier_SO4%cutoff
      this%species_Z => this%fourier_SO4%species_Z
      this%w => this%fourier_SO4%w
      
      this%initialised = .true.

   endsubroutine power_so4_initialise

   subroutine power_so4_finalise(this,error)
      type(power_so4), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      call finalise(this%fourier_SO4,error)

      this%cutoff => null()
      this%z0_ratio => null()
      this%z0 => null()
      this%j_max => null()
      this%Z => null()
      this%cutoff => null()
      this%species_Z => null()
      this%w => null()

      this%initialised = .false.

   endsubroutine power_so4_finalise

   subroutine soap_initialise(this,args_str,error)
      type(soap), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      real(dp) :: alpha_basis, spacing_basis, cutoff_basis, basis_error_exponent
      real(dp), dimension(:,:), allocatable :: covariance_basis, overlap_basis, cholesky_overlap_basis
      integer :: i, j, xml_version

      type(LA_Matrix) :: LA_covariance_basis, LA_overlap_basis
      character(len=STRING_LENGTH) :: species_Z_str
      logical :: has_n_species, has_species_Z, has_central_reference_all_species


      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', PARAM_MANDATORY, this%cutoff, help_string="Cutoff for soap-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Cutoff transition width for soap-type descriptors")
      call param_register(params, 'l_max', PARAM_MANDATORY, this%l_max, help_string="L_max (spherical harmonics basis band limit) for soap-type descriptors")
      call param_register(params, 'n_max', PARAM_MANDATORY, this%n_max, help_string="N_max (number of radial basis functions) for soap-type descriptors")
      call param_register(params, 'atom_sigma', PARAM_MANDATORY, this%atom_sigma, help_string="Width of atomic Gaussians for soap-type descriptors")
      call param_register(params, 'central_weight', '1.0', this%central_weight, help_string="Weight of central atom in environment")
      call param_register(params, 'central_reference_all_species', 'F', this%central_reference_all_species, has_value_target=has_central_reference_all_species, &
           help_string="Place a Gaussian reference for all atom species densities."// &
           "By default (F) only consider when neighbour is the same species as centre")
      call param_register(params, 'average', 'F', this%global, help_string="Whether to calculate averaged SOAP - one descriptor per atoms object. If false (default) atomic SOAP is returned.")
      call param_register(params, 'diagonal_radial', 'F', this%diagonal_radial, help_string="Only return the n1=n2 elements of the power spectrum.")
 
      call param_register(params, 'covariance_sigma0', '0.0', this%covariance_sigma0, help_string="sigma_0 parameter in polynomial covariance function")
      call param_register(params, 'normalise', 'T', this%normalise, help_string="Normalise descriptor so magnitude is 1. In this case the kernel of two equivalent environments is 1.")
      call param_register(params, 'basis_error_exponent', '10.0', basis_error_exponent, help_string="10^(-basis_error_exponent) is the max difference between the target and the expanded function")

      call param_register(params, 'n_Z', '1', this%n_Z, help_string="How many different types of central atoms to consider")
      call param_register(params, 'n_species', '1', this%n_species, has_value_target=has_n_species, help_string="Number of species for the descriptor")
      call param_register(params, 'species_Z', '', species_Z_str, has_value_target=has_species_Z, help_string="Atomic number of species")
      call param_register(params, 'xml_version', '1426512068', xml_version, help_string="Version of GAP the XML potential file was created")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='soap_initialise args_str')) then
         RAISE_ERROR("soap_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      ! backwards compatibility: the default used to be different before this version number
      if( xml_version < 1426512068 ) this%central_reference_all_species = .true.

      allocate(this%species_Z(this%n_species))
      allocate(this%Z(this%n_Z))

      if( has_species_Z .and. .not. has_n_species ) then
         RAISE_ERROR("soap_initialise: is species_Z is present, n_species must be present, too.",error)
      endif

      call initialise(params)
      if( has_n_species ) then
         if(this%n_species == 1) then
            call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         else
            call param_register(params, 'species_Z', '//MANDATORY//', this%species_Z, help_string="Atomic number of species")
         endif
      else
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
      endif

      if( .not. has_central_reference_all_species .and. this%n_species == 1 ) this%central_reference_all_species = .true.

      if( this%n_Z == 1 ) then
         call param_register(params, 'Z', '0', this%Z(1), help_string="Atomic number of central atom, 0 is the wild-card")
      else
         call param_register(params, 'Z', '//MANDATORY//', this%Z, help_string="Atomic numbers to be considered for central atom, must be a list")
      endif
      
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='soap_initialise args_str')) then
         RAISE_ERROR("soap_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%alpha = 0.5_dp / this%atom_sigma**2

      alpha_basis = this%alpha
      cutoff_basis = this%cutoff + this%atom_sigma * sqrt(2.0_dp * basis_error_exponent * log(10.0_dp))
      spacing_basis = cutoff_basis / this%n_max

      allocate(this%r_basis(this%n_max), this%transform_basis(this%n_max,this%n_max), &
         covariance_basis(this%n_max,this%n_max), overlap_basis(this%n_max,this%n_max), this%cholesky_overlap_basis(this%n_max,this%n_max))

      !this%r_basis(this%n_max) = cutoff_basis
      !do i = this%n_max-1, 1, -1
      !   this%r_basis(i)  = this%r_basis(i+1) - spacing_basis
      !enddo

      this%r_basis(1) = 0.0_dp
      do i = 2, this%n_max
         this%r_basis(i)  = this%r_basis(i-1) + spacing_basis
      enddo

      do i = 1, this%n_max
         do j = 1, this%n_max
            covariance_basis(j,i) = exp(-alpha_basis * (this%r_basis(i) - this%r_basis(j))**2)
            !overlap_basis(j,i) = exp(-0.5_dp * alpha_basis* (this%r_basis(i) - this%r_basis(j))**2) * ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) )
            !print*, 'A', exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) )
            !print*, 'B', sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j))
            !print*, 'C', alpha_basis*exp(0.5_dp * alpha_basis * (this%r_basis(i) + this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 )
            !print*, 'D', ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) )
            !overlap_basis(j,i) = exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) ) * &
            !   ( sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j)) + &
            !   alpha_basis*exp(0.5_dp * alpha_basis * (this%r_basis(i) + this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 ) * &
            !   ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) ) )

            overlap_basis(j,i) = ( exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) ) * &
               sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j)) + &
               alpha_basis*exp(-0.5_dp * alpha_basis * (this%r_basis(i) - this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 ) * &
               ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) ) )
         enddo
      enddo

      !overlap_basis = overlap_basis * sqrt(pi / ( 8.0_dp * alpha_basis ) )
      overlap_basis = overlap_basis / sqrt(128.0_dp * alpha_basis**5)

      call initialise(LA_covariance_basis,covariance_basis)
      call initialise(LA_overlap_basis,overlap_basis)
      call LA_Matrix_Factorise(LA_overlap_basis, this%cholesky_overlap_basis)
      do i = 1, this%n_max
         do j = 1, i-1 !i + 1, this%n_max
            this%cholesky_overlap_basis(j,i) = 0.0_dp
         enddo
      enddo

      call Matrix_Solve(LA_covariance_basis,this%cholesky_overlap_basis,this%transform_basis)

      call finalise(LA_covariance_basis)
      call finalise(LA_overlap_basis)

      if(allocated(covariance_basis)) deallocate(covariance_basis)
      if(allocated(overlap_basis)) deallocate(overlap_basis)
      if(allocated(cholesky_overlap_basis)) deallocate(cholesky_overlap_basis)

      this%initialised = .true.

   endsubroutine soap_initialise

   subroutine soap_finalise(this,error)
      type(soap), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%l_max = 0
      this%alpha = 0.0_dp
      this%central_weight = 0.0_dp
      this%central_reference_all_species = .false.
      this%global = .false.
      this%diagonal_radial = .false.
      this%covariance_sigma0 = 0.0_dp
      this%normalise = .true.

      this%n_max = 0
      this%n_Z = 0
      this%n_species = 0

      if(allocated(this%r_basis)) deallocate(this%r_basis)
      if(allocated(this%transform_basis)) deallocate(this%transform_basis)
      if(allocated(this%cholesky_overlap_basis)) deallocate(this%cholesky_overlap_basis)
      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%Z)) deallocate(this%Z)

      this%initialised = .false.

   endsubroutine soap_finalise

   subroutine AN_monomer_initialise(this,args_str,error)
      type(AN_monomer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for AN_monomer-type descriptors")
      call param_register(params, 'atomic_number', '1', this%atomic_number, help_string="Atomic number in AN_monomer-type descriptors")
      call param_register(params, 'N', '4', this%N, help_string="Number of atoms in cluster")
      call param_register(params, 'do_atomic', 'T', this%do_atomic, help_string="Descriptors are cluster based or atom-based")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='AN_monomer_initialise args_str')) then
         RAISE_ERROR("AN_monomer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine AN_monomer_initialise

   subroutine AN_monomer_finalise(this,error)
      type(AN_monomer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%atomic_number = 0
      this%N = 0

      this%do_atomic = .false.
      this%initialised = .false.

   endsubroutine AN_monomer_finalise

   subroutine general_monomer_initialise(this,args_str,error)
      type(general_monomer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_string
      character(len=STRING_LENGTH), dimension(99) :: signature_fields
      integer, optional, intent(out) :: error
      integer :: i,n_atoms,j


      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for general_monomer-type descriptors")
      call param_register(params, 'signature', PARAM_MANDATORY, signature_string, help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'atom_ordercheck', 'true', this%atom_ordercheck, help_string="T: find molecules. F: go by order of atoms")
      call param_register(params, 'strict', 'true', this%strict, help_string="Raise error if not all atoms assigned to monomer")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='general_monomer_initialise args_str')) then
         RAISE_ERROR("general_monomer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call split_string(signature_string,' ','{}',signature_fields(:),n_atoms,matching=.true.)
      allocate(this%signature(n_atoms))

      do i=1,n_atoms
        this%signature(i) = string_to_int(signature_fields(i))
      end do

      call permutation_data_initialise(this%permutation_data,signature_one=this%signature,error=error)

      this%initialised = .true.

   endsubroutine general_monomer_initialise

   subroutine general_monomer_finalise(this,error)
      type(general_monomer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      if(allocated(this%signature)) deallocate(this%signature)

      this%initialised = .false.

   endsubroutine general_monomer_finalise

   subroutine com_dimer_initialise(this,args_str,error)
      type(com_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_one_string, signature_two_string
      character(len=STRING_LENGTH), dimension(99) :: signature_one_fields, signature_two_fields
      integer, optional, intent(out) :: error
      integer :: i, n_atoms_one, n_atoms_two

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff(intermolecular) for com_dimer-type descriptors")
      call param_register(params, 'monomer_one_cutoff', '0.00', this%monomer_one_cutoff, help_string="Cutoff(mono1) for com_dimer-type descriptors")
      call param_register(params, 'monomer_two_cutoff', '0.00', this%monomer_two_cutoff, help_string="Cutoff(mono2) for com_dimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Width of smooth cutoff region for com_dimer-type descriptors")
      call param_register(params, 'atom_ordercheck', 'true', this%atom_ordercheck, help_string="T: find molecules. F: go by order of atoms")
      call param_register(params, 'strict', 'true', this%strict, help_string="Raise error if not all atoms assigned to monomer or if no monomer pairs found")
      call param_register(params, 'mpifind', 'false', this%mpifind, help_string="Use find_monomer_pairs_MPI")
      call param_register(params, 'signature_one', PARAM_MANDATORY, signature_one_string, help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'signature_two', PARAM_MANDATORY, signature_two_string, help_string="Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='com_dimer_initialise args_str')) then
         RAISE_ERROR("com_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%transfer_parameters, args_str, error)

      call split_string(signature_one_string,' ','{}',signature_one_fields(:),n_atoms_one,matching=.true.)
      call split_string(signature_two_string,' ','{}',signature_two_fields(:),n_atoms_two,matching=.true.)
      allocate(this%signature_one(n_atoms_one))
      allocate(this%signature_two(n_atoms_two))

      do i=1,n_atoms_one
        this%signature_one(i) = string_to_int(signature_one_fields(i))
      end do
      do i=1,n_atoms_two
        this%signature_two(i) = string_to_int(signature_two_fields(i))
      end do

      this%monomers_identical=.False.
      if (size(this%signature_one) == size(this%signature_two)) then
         if (all(this%signature_one == this%signature_two)) then
            this%monomers_identical = .True.
         end if
      end if

      this%initialised = .true.
   endsubroutine com_dimer_initialise

   subroutine com_dimer_finalise(this,error)
      type(com_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)
      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_one_cutoff = 0.0_dp
      this%monomer_two_cutoff = 0.0_dp
      this%atom_ordercheck = .true.
      this%use_smooth_cutoff = .false.
      if(allocated(this%signature_one)) deallocate(this%signature_one)
      if(allocated(this%signature_two)) deallocate(this%signature_two)

      this%initialised = .false.

   endsubroutine com_dimer_finalise

   subroutine general_dimer_initialise(this,args_str,error)
      type(general_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_one_string, signature_two_string
      character(len=STRING_LENGTH), dimension(99) :: signature_one_fields, signature_two_fields
      integer, optional, intent(out) :: error
      integer :: i,j, n_atoms_one, n_atoms_two, dimer_size, start, finish, d
      logical, dimension(:,:), allocatable :: intermolecular
      integer, dimension(:), allocatable :: signature

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff(intermolecular) for general_dimer-type descriptors")
      call param_register(params, 'monomer_one_cutoff', '0.00', this%monomer_one_cutoff, help_string="Cutoff(mono1) for general_dimer-type descriptors")
      call param_register(params, 'monomer_two_cutoff', '0.00', this%monomer_two_cutoff, help_string="Cutoff(mono2) for general_dimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Width of smooth cutoff region for general_dimer-type descriptors")
      call param_register(params, 'internal_swaps_only', 'true', this%internal_swaps_only, help_string="F: energies will be symmetrised over swaps of nuclei between monomers")
      call param_register(params, 'atom_ordercheck', 'true', this%atom_ordercheck, help_string="T: find molecules. F: go by order of atoms")
      call param_register(params, 'double_count', 'false', this%double_count, help_string="T: double count when constructing the dimers, for compatibility with water dimer descriptor, default False")
      call param_register(params, 'strict', 'true', this%strict, help_string="Raise error if not all atoms assigned to monomer or if no monomer pairs found")
      call param_register(params, 'use_com', 'false', this%use_com, help_string="Use COM instead of COG")
      call param_register(params, 'mpifind', 'false', this%mpifind, help_string="Use find_monomer_pairs_MPI")
      call param_register(params, 'signature_one', PARAM_MANDATORY, signature_one_string, help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'signature_two', PARAM_MANDATORY, signature_two_string, help_string="Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='general_dimer_initialise args_str')) then
         RAISE_ERROR("general_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%transfer_parameters, args_str, error)

      call split_string(signature_one_string,' ','{}',signature_one_fields(:),n_atoms_one,matching=.true.)
      call split_string(signature_two_string,' ','{}',signature_two_fields(:),n_atoms_two,matching=.true.)
      allocate(this%signature_one(n_atoms_one))
      allocate(this%signature_two(n_atoms_two))

      do i=1,n_atoms_one
        this%signature_one(i) = string_to_int(signature_one_fields(i))
      end do
      do i=1,n_atoms_two
        this%signature_two(i) = string_to_int(signature_two_fields(i))
      end do

      this%monomers_identical=.False.
      if (size(this%signature_one) == size(this%signature_two)) then
         if (all(this%signature_one == this%signature_two)) then
            this%monomers_identical = .True.
         end if
      end if

      call permutation_data_initialise(this%permutation_data,signature_one=this%signature_one,signature_two=this%signature_two,internal_swaps_only=this%internal_swaps_only,error=error)

      dimer_size=n_atoms_one + n_atoms_two
      d=dimer_size*(dimer_size-1)/2

      allocate(signature(dimer_size))
      allocate(intermolecular(dimer_size,dimer_size))
      allocate(this%is_intermolecular(d))
      allocate(this%cutoff_contributor(d))
      allocate(this%component_atoms(d,2))

      signature(1:n_atoms_one) = this%signature_one
      signature(1+n_atoms_one:dimer_size) = this%signature_two
      intermolecular = .false.
      this%cutoff_contributor=.false.

      do i=1,n_atoms_one
        do j=1+n_atoms_one,dimer_size
          intermolecular(i,j)=.true.
        end do
      end do

      start = 0
      finish=dimer_size-1
      do i=1,dimer_size
        do j=1,finish-start
          this%is_intermolecular(start+j) = intermolecular(i,i+j)
          this%component_atoms(start+j,:) = (/ i, i+j /)
        end do
        start = finish
        finish=finish + dimer_size-i-1
      end do


      do i=1,d
        if (this%is_intermolecular(i)) then
          if (.not. signature(this%component_atoms(i,1))==1 ) then
            if (.not. signature(this%component_atoms(i,2))==1 ) then
              this%cutoff_contributor(i)=.true.
            end if
          end if
        end if
      end do

      this%initialised = .true.

      deallocate(signature)
      deallocate(intermolecular)
   endsubroutine general_dimer_initialise

   subroutine general_dimer_finalise(this,error)
      type(general_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)
      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_one_cutoff = 0.0_dp
      this%monomer_two_cutoff = 0.0_dp
      this%atom_ordercheck = .true.
      this%internal_swaps_only = .true.
      this%use_smooth_cutoff = .false.
      if(allocated(this%signature_one)) deallocate(this%signature_one)
      if(allocated(this%signature_two)) deallocate(this%signature_two)
      if(allocated(this%is_intermolecular)) deallocate(this%is_intermolecular)
      if(allocated(this%component_atoms)) deallocate(this%component_atoms)
      if(allocated(this%cutoff_contributor)) deallocate(this%cutoff_contributor)

      this%initialised = .false.

   endsubroutine general_dimer_finalise

   subroutine general_trimer_initialise(this,args_str,error)
      type(general_trimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_one_string, signature_two_string, signature_three_string
      character(len=STRING_LENGTH), dimension(99) :: signature_one_fields, signature_two_fields, signature_three_fields
      integer, optional, intent(out) :: error
      integer :: i,j, n_atoms_one, n_atoms_two, n_atoms_three, trimer_size, start, finish,d
      logical, dimension(:,:), allocatable :: intermolecular
      integer, dimension(:), allocatable :: signature

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff(intermolecular) for general_trimer-type descriptors")
      call param_register(params, 'monomer_one_cutoff', '0.00', this%monomer_one_cutoff, help_string="Cutoff(mono1) for general_trimer-type descriptors")
      call param_register(params, 'monomer_two_cutoff', '0.00', this%monomer_two_cutoff, help_string="Cutoff(mono2) for general_trimer-type descriptors")
      call param_register(params, 'monomer_three_cutoff', '0.00', this%monomer_three_cutoff, help_string="Cutoff(mono3) for general_trimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Width of smooth cutoff region for general_trimer-type descriptors")
      call param_register(params, 'internal_swaps_only', 'true', this%internal_swaps_only, help_string="F: energies will be symmetrised over swaps of nuclei between monomers")
      call param_register(params, 'atom_ordercheck', 'true', this%atom_ordercheck, help_string="T: find molecules. F: go by order of atoms")
      call param_register(params, 'strict', 'true', this%strict, help_string="Raise error if not all atoms assigned to monomer or if no monomer pairs found")
      call param_register(params, 'use_com', 'false', this%use_com, help_string="Use COM instead of COG")
      call param_register(params, 'mpifind', 'false', this%mpifind, help_string="Use find_monomer_triplets_MPI")
      call param_register(params, 'signature_one', PARAM_MANDATORY, signature_one_string, help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'signature_two', PARAM_MANDATORY, signature_two_string, help_string="Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'signature_three', PARAM_MANDATORY, signature_three_string, help_string="Atomic numbers of monomer three, format {Z1 Z2 Z3 ...}")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='general_trimer_initialise args_str')) then
         RAISE_ERROR("general_trimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call split_string(signature_one_string,' ','{}',signature_one_fields(:),n_atoms_one,matching=.true.)
      call split_string(signature_two_string,' ','{}',signature_two_fields(:),n_atoms_two,matching=.true.)
      call split_string(signature_three_string,' ','{}',signature_three_fields(:),n_atoms_three,matching=.true.)
      allocate(this%signature_one(n_atoms_one))
      allocate(this%signature_two(n_atoms_two))
      allocate(this%signature_three(n_atoms_three))

      do i=1,n_atoms_one
        this%signature_one(i) = string_to_int(signature_one_fields(i))
      end do
      do i=1,n_atoms_two
        this%signature_two(i) = string_to_int(signature_two_fields(i))
      end do
      do i=1,n_atoms_three
        this%signature_three(i) = string_to_int(signature_three_fields(i))
      end do

      this%one_two_identical = .false.
      this%one_three_identical = .false.
      this%two_three_identical = .false.

      if (size(this%signature_one) == size(this%signature_two)) then
         if (all(this%signature_one == this%signature_two)) then
            this%one_two_identical = .True.
         end if
      end if

      if (size(this%signature_one) == size(this%signature_three)) then
         if (all(this%signature_one == this%signature_three)) then
            this%one_three_identical = .True.
         end if
      end if

      if (size(this%signature_two) == size(this%signature_three)) then
         if (all(this%signature_two == this%signature_three)) then
            this%two_three_identical = .True.
         end if
      end if

      call permutation_data_initialise(this%permutation_data,signature_one=this%signature_one,signature_two=this%signature_two,signature_three=this%signature_three,internal_swaps_only=this%internal_swaps_only,error=error)

      trimer_size=n_atoms_one + n_atoms_two + n_atoms_three
      d=trimer_size*(trimer_size-1)/2

      allocate(signature(trimer_size))
      allocate(intermolecular(trimer_size,trimer_size))
      allocate(this%is_intermolecular(d))
      allocate(this%cutoff_contributor(d))
      allocate(this%component_atoms(d,2))

      signature(1:n_atoms_one) = this%signature_one
      signature(1+n_atoms_one:n_atoms_one+n_atoms_two) = this%signature_two
      signature(1+n_atoms_one+n_atoms_two:trimer_size) = this%signature_three
      intermolecular = .false.
      this%cutoff_contributor=.false.

      do i=1,n_atoms_one
        do j=1+n_atoms_one,trimer_size
          intermolecular(i,j)=.true.
          intermolecular(j,i)=.true.
        end do
      end do
      do i=1+n_atoms_one,n_atoms_one+n_atoms_two
        do j=1+n_atoms_one+n_atoms_two,trimer_size
          intermolecular(i,j)=.true.
          intermolecular(j,i)=.true.
        end do
      end do

      start = 0
      finish=trimer_size-1
      do i=1,trimer_size
        do j=1,finish-start
          this%is_intermolecular(start+j) = intermolecular(i,i+j)
          this%component_atoms(start+j,:) = (/ i, i+j /)
        end do
        start = finish
        finish=finish + trimer_size-i-1
      end do

      do i=1,d
        if (this%is_intermolecular(i)) then
          if (.not. signature(this%component_atoms(i,1))==1 ) then
            if (.not. signature(this%component_atoms(i,2))==1 ) then
              this%cutoff_contributor(i)=.true.
            end if
          end if
        end if
      end do

      this%initialised = .true.
      deallocate(signature)
      deallocate(intermolecular)
   endsubroutine general_trimer_initialise

   subroutine general_trimer_finalise(this,error)
      type(general_trimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)
      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_one_cutoff = 0.0_dp
      this%monomer_two_cutoff = 0.0_dp
      this%monomer_three_cutoff = 0.0_dp
      this%atom_ordercheck = .true.
      this%internal_swaps_only = .true.
      this%use_smooth_cutoff = .false.
      if(allocated(this%signature_one)) deallocate(this%signature_one)
      if(allocated(this%signature_two)) deallocate(this%signature_two)
      if(allocated(this%signature_three)) deallocate(this%signature_three)

      this%initialised = .false.

   endsubroutine general_trimer_finalise

   subroutine rdf_initialise(this,args_str,error)
      type(rdf), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: i
      real(dp) :: r_min, r_max
      logical :: has_r_max, has_w_gauss

      INIT_ERROR(error)

      call finalise(this)
      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for rdf-type descriptors")
      call param_register(params, 'transition_width', '0.20', this%transition_width, help_string="Width of transition region from 1 to 0")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'r_min', '0.0', r_min, help_string="Atomic number of central atom")
      call param_register(params, 'r_max', '0.0', r_max, has_value_target = has_r_max, help_string="Atomic number of central atom")
      call param_register(params, 'n_gauss', '10', this%n_gauss, help_string="Atomic number of central atom")
      call param_register(params, 'w_gauss', '0.0', this%w_gauss, has_value_target = has_w_gauss, help_string="Atomic number of central atom")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='rdf_initialise args_str')) then
         RAISE_ERROR("rdf_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%r_gauss(this%n_gauss))
      if(.not. has_w_gauss) this%w_gauss = this%cutoff / this%n_gauss * 2.0_dp
      if(.not. has_r_max) r_max = this%cutoff - this%w_gauss / 2.0_dp
      this%r_gauss = real( (/(i,i=1,this%n_gauss)/), kind=dp ) / real(this%n_gauss,kind=dp) * (r_max - r_min) + r_min

      this%initialised = .true.

   endsubroutine rdf_initialise

   subroutine rdf_finalise(this,error)
      type(rdf), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%transition_width = 0.0_dp
      this%Z = 0
      this%n_gauss = 0
      if( allocated(this%r_gauss) ) deallocate(this%r_gauss)

      this%initialised = .false.

   endsubroutine rdf_finalise

   subroutine as_distance_2b_initialise(this,args_str,error)
      type(as_distance_2b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Lower cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'max_cutoff', PARAM_MANDATORY, this%max_cutoff, help_string="Higher cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'as_cutoff', PARAM_MANDATORY, this%as_cutoff, help_string="Cutoff of asymmetricity")
      call param_register(params, 'overlap_alpha', '0.50', this%as_cutoff, help_string="Cutoff of asymmetricity")
      call param_register(params, 'min_transition_width', '0.50', this%min_transition_width, help_string="Transition width of lower cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'max_transition_width', '0.50', this%max_transition_width, help_string="Transition width of higher cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'as_transition_width', '0.10', this%as_transition_width, help_string="Transition width of asymmetricity cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'coordination_cutoff', PARAM_MANDATORY, this%coordination_cutoff, help_string="Cutoff for coordination function in as_distance_2b-type descriptors")
      call param_register(params, 'coordination_transition_width', '0.50', this%coordination_transition_width, help_string="Transition width for as_distance_2b-type descriptors")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atom type #1 in bond")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atom type #2 in bond")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='as_distance_2b_initialise args_str')) then
         RAISE_ERROR("as_distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine as_distance_2b_initialise

   subroutine as_distance_2b_finalise(this,error)
      type(as_distance_2b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%min_cutoff = 0.0_dp
      this%max_cutoff = 0.0_dp
      this%as_cutoff = 0.0_dp
      this%overlap_alpha = 0.0_dp
      this%min_transition_width = 0.0_dp
      this%max_transition_width = 0.0_dp
      this%as_transition_width = 0.0_dp
      this%coordination_cutoff = 0.0_dp
      this%coordination_transition_width = 0.0_dp
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine as_distance_2b_finalise

   subroutine molecule_lo_d_initialise(this,args_str,error)
      type(molecule_lo_d), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_string, atoms_template_string, symmetry_string, symmetry_property_name, append_file, append_string
      character(len=STRING_LENGTH), dimension(99) :: signature_fields, symmetry_rows, row_fields, append_rows,template_rows
      integer, optional, intent(out) :: error
      integer :: i,n_atoms,j,n_symm_rows, current_depth, start, finish, i_component, atom_j, N_atom_pairs, atom_k,n_append_rows, old_size, n_perms,n_template_rows
      integer, dimension(:,:), allocatable :: equivalents_input, bonds_to_append, tmp_permutations
      logical :: signature_given, symmetries_given, add_bond, j_k_present, k_j_present
      integer, dimension(:,:), pointer :: symm_2d
      integer, dimension(:), pointer :: symm_1d

      type(Table) :: atom_a, atom_b
      type(CInOutput) :: tempatoms
      type(inoutput) :: tempfile
      type(inoutput) :: symmetry_inout
      type(inoutput) :: append_inout
      type(Connection) :: at_connect

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for molecule_lo_d-type descriptors")
      call param_register(params, 'atoms_template_string', PARAM_MANDATORY, atoms_template_string , help_string="Atoms object which serves as a template - written to a string")
      call param_register(params, 'neighbour_graph_depth', '2', this%neighbour_graph_depth, help_string="Ignore distances between atoms separated by more than this number of bonds")
      call param_register(params, 'signature', '', signature_string, help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'symmetry_property_name', 'symm', symmetry_property_name, help_string="Integer arrays specifying symmetries - see header of make_permutations_v2.f95 for format")
      call param_register(params, 'append_file', '', append_file, help_string="Pairs of atoms for which we want the distance to be additionally included in the descriptor")
      call param_register(params, 'atom_ordercheck', 'T', this%atom_ordercheck, help_string= &
                                                                    "T: basic check that atoms in same order as in template F: assume all xyz frames have atoms in same order")
      call param_register(params, 'desctype', '0', this%desctype, help_string="0: distance matrix, 1: inverse distance matrix, 2: Coulomb matrix, 3: exponential")
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='molecule_lo_d_initialise args_str')) then
         RAISE_ERROR("molecule_lo_d_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif

      call finalise(params)

      do i=1,len_trim(atoms_template_string)
        if(atoms_template_string(i:i)=='%') then
          atoms_template_string(i:i)=' '
        end if
      end do

      ! read in the atoms object in the sample geometry file
      call initialise(tempfile, filename="temp.xyz",action=OUTPUT)
      call split_string(atoms_template_string,';','{}',template_rows(:),n_template_rows,matching=.true.)
      do i=1,n_template_rows
        call print(template_rows(i),file=tempfile)
      end do
      call finalise(tempfile)
      call initialise(tempatoms,"temp.xyz")
      call read(this%template_atoms, tempatoms, error=error)

      this%n_atoms = this%template_atoms%N
      ! make a table of bonds - this copied from topology module private function create_bond_list
      call calc_connect(at_connect,this%template_atoms,error=error)
      if (this%neighbour_graph_depth > 0) then
        do i=1,this%template_atoms%N

           call initialise(atom_a,4,0,0,0,0)
           call append(atom_a,(/i,0,0,0/))
           call bfs_step(this%template_atoms,atom_a,atom_b,nneighb_only=.false.,min_images_only=.true.,alt_connect=at_connect)

           do j = 1,atom_b%N
              atom_j = atom_b%int(1,j)
              if (atom_j.gt.i) then
                 add_bond = .true.

                 if (add_bond) then
                    call append(this%bonds,(/i,atom_j/))
                 endif

              else
!                 call print('not added '//i//' -- '//atom_j)
              endif
           enddo
           call finalise(atom_a)
           call finalise(atom_b)
        enddo
        ! add atom pairs separated by up to neighbour_graph_depth bonds
        current_depth = 1
        this%atom_pairs=this%bonds
        do while (current_depth < this%neighbour_graph_depth)
           current_depth = current_depth + 1
           call bond_list_next_layer(this%bonds,this%atom_pairs)
        end do
      endif

     ! append any manually specified bonds, if not already present
      if (.not. append_file .eq. '') then
         call initialise(append_inout,trim(append_file))
         read(append_inout%unit,'(a)') append_string
         call split_string(append_string,',','{}',append_rows(:),n_append_rows,matching=.true.)
         
         do i = 1, n_append_rows
            call split_string(append_rows(i),' ','{}',row_fields(:),n_atoms,matching=.true.)
            if ( n_atoms .ne. 2) then
              RAISE_ERROR("append_file incorrectly formatted, expected atoms in pairs",error)
            end if
            if (.not. allocated(bonds_to_append)) allocate(bonds_to_append(n_append_rows,n_atoms))
            do j=1,n_atoms
              bonds_to_append(i,j) = string_to_int(row_fields(j))
            end do
         end do


         do i=1,size(bonds_to_append,1)
           atom_j = bonds_to_append(i,1)
           atom_k = bonds_to_append(i,2)
           if ( this%neighbour_graph_depth > 0) then
             N_atom_pairs = this%atom_pairs%N
             j_k_present = any(this%atom_pairs%int(1,:N_atom_pairs) .eq. atom_j .and. (this%atom_pairs%int(2,:N_atom_pairs) .eq. atom_k))
             k_j_present = any(this%atom_pairs%int(1,:N_atom_pairs) .eq. atom_k .and. (this%atom_pairs%int(2,:N_atom_pairs) .eq. atom_j))
           else
             j_k_present = .false.
             k_j_present = .false.
           end if
           if (.not. j_k_present .and. .not. k_j_present) then
             call append(this%atom_pairs,(/atom_j,atom_k/))
           end if
         end do
      end if

      !call print('table of atom pairs included in descriptor')
      !call print(this%atom_pairs)
      signature_given = .False.
      symmetries_given = .False.
      call permutation_data_finalise(this%permutation_data)
      ! parse signature, if given
      if (.not. signature_string .eq. '') then
         signature_given=.True.
         call split_string(signature_string,' ','{}',signature_fields(:),n_atoms,matching=.true.)
         if (.not. n_atoms .eq. this%n_atoms) then
           RAISE_ERROR('signature does not have correct number of entries', error)
         end if
         allocate(this%signature(this%n_atoms))

         do i=1,this%n_atoms
           this%signature(i) = string_to_int(signature_fields(i))
         end do
         call permutation_data_initialise(this%permutation_data,signature_one=this%signature,error=error)
      end if

      ! parse symmetries, if given
      symmetries_given = has_property(this%template_atoms,symmetry_property_name)
      if (symmetries_given) then
        if (.not. assign_pointer(this%template_atoms, symmetry_property_name, symm_1d)) then
           if (.not. assign_pointer(this%template_atoms, symmetry_property_name, symm_2d)) then
             RAISE_ERROR('IPModel_Coulomb_Calc failed to assign pointer to "'//trim(symmetry_property_name)//'" property', error)
           else
             if (.not. allocated(equivalents_input)) allocate(equivalents_input(size(symm_2d,1),size(symm_2d,2)))
             equivalents_input=symm_2d
           end if
        else
          if (.not. allocated(equivalents_input)) allocate(equivalents_input(1,size(symm_1d)))
             equivalents_input(1,:)=symm_1d
        endif
        call permutation_data_initialise(this%permutation_data,equivalents_input=equivalents_input,error=error)
      end if

      ! If no signature is given and no symmetries are manually specified we just make a dummy signature to initialise the permutation data
      if (.not. (signature_given) .and. .not. (symmetries_given)) then
         allocate(this%signature(this%n_atoms))
         do i=1,this%n_atoms
           this%signature(i) = i
         end do
         call permutation_data_initialise(this%permutation_data,signature_one=this%signature,error=error)
      end if

      ! make a mapping from i,j pairs to dist_vec components, and compile the list of these which actually make up the descriptor

      this%max_dimension = this%n_atoms * (this%n_atoms -1 ) / 2
      N_atom_pairs = this%atom_pairs%N
      allocate(this%component_atoms(this%max_dimension,2))
      allocate(this%included_components(N_atom_pairs))
      this%included_components=0

      i_component=1
      start = 0
      finish=this%n_atoms-1
      do i=1,this%n_atoms
        do j=1,finish-start
          this%component_atoms(start+j,:) = (/ i, i+j /)
          if (any(this%atom_pairs%int(1,:n_atom_pairs) .eq. i .and. this%atom_pairs%int(2,:n_atom_pairs) .eq. i+j) .or. &
              any(this%atom_pairs%int(2,:n_atom_pairs) .eq. i .and. this%atom_pairs%int(1,:n_atom_pairs) .eq. i+j))  then
               this%included_components(i_component) = start+j
               i_component = i_component + 1
          end if
        end do
        start = finish
        finish=finish + this%n_atoms-i-1
      end do

      if (any(this%included_components .eq. 0)) then
        RAISE_ERROR('molecule_lo_d_initialise : something went wrong picking out the correct interatomic distances',error)
      end if

      this%initialised = .true.

   endsubroutine molecule_lo_d_initialise

   subroutine molecule_lo_d_finalise(this,error)
      type(molecule_lo_d), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      if (allocated(this%signature)) deallocate(this%signature)
      if (allocated(this%component_atoms)) deallocate(this%component_atoms)
      if (allocated(this%included_components)) deallocate(this%included_components)



      this%initialised = .false.

   endsubroutine molecule_lo_d_finalise

   subroutine alex_initialise(this,args_str,error)
      type(alex), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for alex-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'power_min', '5', this%power_min, help_string="Minimum power of radial basis for the descriptor")
      call param_register(params, 'power_max', '10', this%power_max, help_string="Maximum power of the radial basis for the descriptor")
      call param_register(params, 'n_species', '1', this%n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='alex_initialise args_str')) then
         RAISE_ERROR("alex_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(this%n_species))

      call initialise(params)
      if( this%n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='alex_initialise args_str')) then
         RAISE_ERROR("alex_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine alex_initialise

   subroutine alex_finalise(this,error)
      type(alex), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp

      if(allocated(this%species_Z)) deallocate(this%species_Z)

      this%initialised = .false.

   endsubroutine alex_finalise

   subroutine distance_Nb_initialise(this,args_str,error)
      type(distance_Nb), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(len=STRING_LENGTH) :: default_Z = ""
      integer :: i, j, k, i_p
      integer :: nEdges, nConnectivities, nMonomerConnectivities
      integer, dimension(:), allocatable :: n_permutations, connectivityList
      integer, dimension(:,:), allocatable :: atom_permutations, distance_matrix_index, edges

      logical, dimension(:,:,:), allocatable :: allConnectivities


      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', PARAM_MANDATORY, this%cutoff, help_string="Cutoff for distance_Nb-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.5', this%cutoff_transition_width, help_string="Transition width of cutoff for distance_Nb-type descriptors")
      call param_register(params, 'order', PARAM_MANDATORY, this%order, help_string="Many-body order, in terms of number of neighbours")
      call param_register(params, 'compact_clusters', "F", this%compact_clusters, help_string="If true, generate clusters where the atoms have at least one connection to the central atom. If false, only clusters where all atoms are connected are generated.")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='distance_Nb_initialise args_str')) then
         RAISE_ERROR("distance_Nb_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if( this%order < 1 ) then
         RAISE_ERROR("distance_Nb_initialise: order must be greater than 0",error)
      endif

      allocate(this%Z(this%order))
      default_Z = ""
      do i = 1, this%order
         default_Z = trim(default_Z) // " 0"
      enddo

      call initialise(params)
      if( this%order == 1 ) then
         call param_register(params, 'Z', trim(default_Z), this%Z(1), help_string="Atomic type of neighbours")
      else
         call param_register(params, 'Z', trim(default_Z), this%Z, help_string="Atomic type of neighbours")
      endif
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='distance_Nb_initialise args_str')) then
         RAISE_ERROR("distance_Nb_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call sort_array(this%Z)
      call distance_Nb_n_permutations(this%Z, n_permutations)
      this%n_permutations = product(factorial_int(n_permutations))

      allocate(atom_permutations(this%order,this%n_permutations))
      call distance_Nb_permutations(n_permutations,atom_permutations)

      allocate(distance_matrix_index(this%order,this%order))
      allocate(this%permutations( max(1,(this%order - 1) * this%order / 2), this%n_permutations))

      if( this%order == 1 ) then
         this%permutations = 1
      else
         k = 0
         do i = 1, this%order
            do j = i+1, this%order
               k = k + 1
               distance_matrix_index(j,i) = k
               distance_matrix_index(i,j) = k
            enddo
         enddo

         do i_p = 1, this%n_permutations
            k = 0
            do i = 1, this%order
               do j = i+1, this%order
                  k = k + 1
                  this%permutations(k,i_p) = distance_matrix_index(atom_permutations(j,i_p), atom_permutations(i,i_p))
               enddo
            enddo
         enddo
      endif

      nEdges = this%order * (this%order - 1) / 2
      allocate( edges(2,nEdges))

      k = 0
      do i = 1, this%order
         do j = i+1, this%order
            k = k + 1
            edges(:,k) = (/i,j/)
         enddo
      enddo

      nConnectivities = 2**nEdges

      allocate(allConnectivities(this%order,this%order,nConnectivities))
      allocate(connectivityList(nEdges))

      nMonomerConnectivities = 0
      do i = 1, nConnectivities
         call integerDigits(i-1,2,connectivityList)
         allConnectivities(:,:,i) = .false.
         do j = 1, nEdges
            allConnectivities(edges(1,j),edges(2,j),i) = ( connectivityList(j) == 1 )
            allConnectivities(edges(2,j),edges(1,j),i) = ( connectivityList(j) == 1 )
         enddo

         if( graphIsConnected( allConnectivities(:,:,i) ) ) nMonomerConnectivities = nMonomerConnectivities + 1
      enddo

      allocate(this%monomerConnectivities(this%order,this%order,nMonomerConnectivities))
      j = 0
      do i = 1, nConnectivities
         if( graphIsConnected( allConnectivities(:,:,i) ) ) then
            j = j + 1
            this%monomerConnectivities(:,:,j) = allConnectivities(:,:,i)
         endif
      enddo

      if(allocated(n_permutations)) deallocate(n_permutations)
      if(allocated(atom_permutations)) deallocate(atom_permutations)
      if(allocated(distance_matrix_index)) deallocate(distance_matrix_index)
      if(allocated(edges)) deallocate(edges)
      if(allocated(allConnectivities)) deallocate(allConnectivities)
      if(allocated(connectivityList)) deallocate(connectivityList)
      this%initialised = .true.

   endsubroutine distance_Nb_initialise

   subroutine distance_Nb_finalise(this,error)
      type(distance_Nb), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.5_dp
      this%order = 0
      this%n_permutations = 0
      this%compact_clusters = .false.
      if(allocated(this%Z)) deallocate(this%Z)
      if(allocated(this%permutations)) deallocate(this%permutations)
      if(allocated(this%monomerConnectivities)) deallocate(this%monomerConnectivities)

      this%initialised = .false.

   endsubroutine distance_Nb_finalise

   subroutine distance_Nb_n_permutations(Z,n_permutations,error)
      integer, dimension(:), intent(in) :: Z
      integer, dimension(:), allocatable :: n_permutations
      integer, optional, intent(out) :: error

      integer :: i
      integer, dimension(:), allocatable :: uniq_Z

      INIT_ERROR(error)

      call uniq(Z,uniq_Z)
      call reallocate(n_permutations,size(uniq_Z))
     
      do i = 1, size(uniq_Z)
         n_permutations(i) = count( uniq_Z(i) == Z )
      enddo
      
      if(allocated(uniq_Z)) deallocate(uniq_Z)

   endsubroutine distance_Nb_n_permutations

   recursive subroutine distance_Nb_permutations(n_permutations,permutations)
      integer, dimension(:), intent(in) :: n_permutations
      integer, dimension(sum(n_permutations),product(factorial_int(n_permutations))), intent(inout) :: permutations
 
      integer, dimension(:), allocatable, save :: current_permutation
      integer :: i, j, n_lo, n_hi
      integer, save :: recursion_level = 0, i_current_permutation = 0
 
      recursion_level = recursion_level + 1
 

      if( recursion_level == 1 ) then
         i_current_permutation = 0
         allocate(current_permutation(sum(n_permutations)))
         current_permutation = 0
      endif
 

      do i = 1, size(n_permutations)
         if( i == 1 ) then
            n_lo = 1
         else
            n_lo = sum(n_permutations(1:i-1)) + 1
         endif
         n_hi = sum(n_permutations(1:i))
         do j = n_lo, n_hi
            if( i_current_permutation < size(permutations,2) ) then
               if( .not. any(j==current_permutation) .and. recursion_level >= n_lo .and. recursion_level <= n_hi ) then
 
                  current_permutation(recursion_level) = j
                  if( recursion_level == sum(n_permutations) ) then
                     i_current_permutation = i_current_permutation + 1
                     permutations(:,i_current_permutation) = current_permutation
                  else
                     call distance_Nb_permutations(n_permutations,permutations)
                  endif
               endif
            endif
         enddo
      enddo
 
      current_permutation(recursion_level) = 0
 
      recursion_level = recursion_level - 1
 
      if( recursion_level == 0 ) then
         deallocate(current_permutation)
      endif
 
   endsubroutine distance_Nb_permutations

   subroutine descriptor_str_add_species(this,species,descriptor_str,error)
      character(len=*), intent(in) :: this
      integer, dimension(:), intent(in) :: species
      character(len=STRING_LENGTH), dimension(:), allocatable, intent(out) :: descriptor_str
      integer, optional, intent(out) :: error

      integer :: my_descriptor_type, i, j, k, l, n_species, order, n
      integer, dimension(:,:), allocatable :: ZN
      real(dp), dimension(:), allocatable :: w
      type(Dictionary) :: params

      INIT_ERROR(error)

      if(allocated(descriptor_str)) deallocate(descriptor_str)

      my_descriptor_type = get_descriptor_type(this,error)
      n_species = size(species)

      select case(my_descriptor_type)
      case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_COSNX,DT_POWER_SO3,DT_POWER_SO4)
         allocate(w(n_species))
         allocate(descriptor_str(n_species))

         if( n_species == 1 ) then
            w = 1.0_dp
         else
            w = real( (/ (i, i=0, n_species-1) /), kind=dp ) / (n_species-1) * 0.5_dp + 0.5_dp
         endif

         do i = 1, n_species
            descriptor_str(i) = trim(this)//" n_species="//n_species//" Z="//species(i)//" species_Z={"//species//"} w={"//w//"}"
         enddo

         deallocate(w)
      case(DT_SOAP)
         allocate(descriptor_str(n_species))
         do i = 1, n_species
            descriptor_str(i) = trim(this)//" n_species="//n_species//" Z="//species(i)//" species_Z={"//species//"}"
         enddo
      case(DT_DISTANCE_2B,DT_CO_DISTANCE_2B,DT_AS_DISTANCE_2B)
         allocate(descriptor_str(n_species * (n_species+1) / 2))

         l = 0
         do i = 1, n_species
            do j = i, n_species
               l = l + 1
               descriptor_str(l) = trim(this)//" Z1="//species(i)//" Z2="//species(j)
            enddo
         enddo

      case(DT_COORDINATION,DT_RDF)
         allocate(descriptor_str(n_species))
         do i = 1, n_species
            descriptor_str(i) = trim(this)//" Z="//species(i)
         enddo
      case(DT_ANGLE_3B,DT_CO_ANGLE_3B)
         allocate(descriptor_str(n_species * n_species * (n_species+1) / 2))
         l = 0
         do i = 1, n_species
            do j = 1, n_species
               do k = j, n_species
                  l = l + 1
                  descriptor_str(l) = trim(this)//" Z="//species(i)//" Z1="//species(j)//" Z2="//species(k)
               enddo
            enddo
         enddo
      case(DT_GENERAL_MONOMER,DT_GENERAL_DIMER,DT_WATER_MONOMER,DT_WATER_DIMER,DT_A2_DIMER,DT_AB_DIMER,DT_TRIHIS,DT_BOND_REAL_SPACE,DT_ATOM_REAL_SPACE,DT_AN_MONOMER)
         allocate(descriptor_str(1))
         descriptor_str(1) = trim(this)
      case(DT_DISTANCE_NB)
         call initialise(params)
         call param_register(params, 'order', PARAM_MANDATORY, order, help_string="Many-body order, in terms of number of neighbours")
         if (.not. param_read_line(params, this, ignore_unknown=.true.,task='descriptor_str_add_species this')) then
            RAISE_ERROR("descriptor_str_add_species failed to parse descriptor string='"//trim(this)//"'", error)
         endif
         call finalise(params)

         n = 1
         do i = 1, order
            n = n * ( n_species + i - 1 ) / i ! avoids double counting
         enddo

         allocate(ZN(order,n),descriptor_str(n))

         call descriptor_str_add_species_distance_Nb(ZN,species,order)

         do i = 1, n
            descriptor_str(i) = trim(this)//" Z={"//ZN(:,i)//"}"
         enddo
         deallocate(ZN)

      case default
         RAISE_ERROR("descriptor_str_add_species: unknown descriptor type "//my_descriptor_type,error)
      endselect

   endsubroutine descriptor_str_add_species

   recursive subroutine descriptor_str_add_species_distance_Nb(ZN,species,order)
      integer, dimension(:,:), intent(inout) :: ZN
      integer, dimension(:), intent(in) :: species
      integer, intent(in) :: order
 
      integer :: i_species, n_species
      integer, save :: current_descriptor, current_order = 0
      integer, dimension(:), allocatable, save :: ZN_current
 
      n_species = size(species)
 
      if( current_order == 0 ) then                           ! first run, outermost order.
         current_descriptor = 0                               ! keeps track of descriptor
         current_order = 1                                    ! keeps track of order
         allocate(ZN_current(order))                          ! builds/updates atomic numbers gradually for each descriptor
      endif
 
      do i_species = 1, n_species
         if( current_order > 1 ) then                             ! no special atom, all atoms equivalent
            if( species(i_species) < ZN_current(current_order-1) ) cycle   ! avoids double-counting of neighbours
         endif
 
         ZN_current(current_order) = species(i_species)
         if( current_order < order ) then                                 ! calls recursively until we reach the last order
            current_order = current_order + 1
            call descriptor_str_add_species_distance_Nb(ZN,species,order)
         else                                                             ! when we reached the last order, fill the atomic numbers in the loop
            current_descriptor = current_descriptor + 1                   ! and add them to the output array
            ZN(:,current_descriptor) = ZN_current
         endif
      enddo
 
      current_order = current_order - 1                                   ! when the loop finished, step one level down
 
      if( current_order == 0 ) deallocate(ZN_current)                     ! when we reach zero, we finished.
 
   endsubroutine descriptor_str_add_species_distance_Nb

   subroutine descriptor_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(cutoff(this) > at%cutoff) then
         RAISE_ERROR("descriptor_calc: descriptor cutoff ("//cutoff(this)//") larger than atoms cutoff("//at%cutoff//")",error)
      endif

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call calc(this%descriptor_bispectrum_SO4,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_BISPECTRUM_SO3)
            call calc(this%descriptor_bispectrum_SO3,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error=error)
         case(DT_BEHLER)
            call calc(this%descriptor_behler,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error=error)
         case(DT_DISTANCE_2b)
            call calc(this%descriptor_distance_2b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_COORDINATION)
            call calc(this%descriptor_coordination,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_ANGLE_3B)
            call calc(this%descriptor_angle_3b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_CO_ANGLE_3B)
            call calc(this%descriptor_co_angle_3b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_CO_DISTANCE_2b)
            call calc(this%descriptor_co_distance_2b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_COSNX)
            call calc(this%descriptor_cosnx,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_TRIHIS)
            call calc(this%descriptor_trihis,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_WATER_MONOMER)
            call calc(this%descriptor_water_monomer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_WATER_DIMER)
            call calc(this%descriptor_water_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_A2_DIMER)
            call calc(this%descriptor_A2_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_AB_DIMER)
            call calc(this%descriptor_AB_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_BOND_REAL_SPACE)
            call calc(this%descriptor_bond_real_space,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_ATOM_REAL_SPACE)
            call calc(this%descriptor_atom_real_space,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_POWER_SO3)
            call calc(this%descriptor_power_so3,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_POWER_SO4)
            call calc(this%descriptor_power_so4,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_SOAP)
            call calc(this%descriptor_soap,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_AN_MONOMER)
            call calc(this%descriptor_AN_monomer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_GENERAL_MONOMER)
            call calc(this%descriptor_general_monomer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_GENERAL_DIMER)
            call calc(this%descriptor_general_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_GENERAL_TRIMER)
            call calc(this%descriptor_general_trimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_RDF)
            call calc(this%descriptor_rdf,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_AS_DISTANCE_2b)
            call calc(this%descriptor_as_distance_2b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_MOLECULE_LO_D)
            call calc(this%descriptor_molecule_lo_d,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_ALEX)
            call calc(this%descriptor_alex,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_COM_DIMER)
            call calc(this%descriptor_com_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_DISTANCE_Nb)
            call calc(this%descriptor_distance_Nb,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case default
            RAISE_ERROR("descriptor_calc: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_calc

   subroutine descriptor_calc_array(this,at,descriptor_out,covariance_cutoff,grad_descriptor_out,grad_descriptor_index,grad_descriptor_pos,grad_covariance_cutoff,args_str,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      real(dp), dimension(:,:), intent(out), optional :: descriptor_out
      real(dp), dimension(:), intent(out), optional :: covariance_cutoff
      real(dp), dimension(:,:,:), intent(out), optional :: grad_descriptor_out
      integer, dimension(:,:), intent(out), optional :: grad_descriptor_index
      real(dp), dimension(:,:), intent(out), optional :: grad_descriptor_pos
      real(dp), dimension(:,:), intent(out), optional :: grad_covariance_cutoff
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(descriptor_data) :: my_descriptor_data
      integer :: i, n, i_d, n_descriptors, n_cross
      logical :: do_grad_descriptor, do_descriptor

      INIT_ERROR(error)

      do_descriptor = present(descriptor_out)
      do_grad_descriptor = present(grad_descriptor_out) .or. present(grad_descriptor_index) .or. present(grad_descriptor_pos)

      call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)

      call calc(this,at,my_descriptor_data,do_descriptor=do_descriptor,do_grad_descriptor=do_grad_descriptor,args_str=args_str,error=error)

      if(present(descriptor_out)) &
         call check_size('descriptor_out',descriptor_out, (/descriptor_dimensions(this),n_descriptors/),'descriptor_calc_array',error)

      if(present(covariance_cutoff)) &
         call check_size('covariance_cutoff',covariance_cutoff,(/n_descriptors/),'descriptor_calc_array',error)

      if(present(grad_descriptor_out)) &
         call check_size('grad_descriptor_out',grad_descriptor_out,(/descriptor_dimensions(this),3,n_cross/),'descriptor_calc_array',error)

      if(present(grad_descriptor_index)) &
         call check_size('grad_descriptor_index',grad_descriptor_index,(/2,n_cross/),'descriptor_calc_array',error)

      if(present(grad_descriptor_pos)) &
         call check_size('grad_descriptor_pos',grad_descriptor_pos,(/3,n_cross/),'descriptor_calc_array',error)

      if(present(grad_covariance_cutoff)) &
         call check_size('grad_covariance_cutoff',grad_covariance_cutoff,(/3,n_cross/),'descriptor_calc_array',error)

      if(do_descriptor) then
         do i = 1, n_descriptors
            descriptor_out(:,i) = my_descriptor_data%x(i)%data
            if(present(covariance_cutoff)) covariance_cutoff(i) = my_descriptor_data%x(i)%covariance_cutoff
         enddo
      endif

      if(do_grad_descriptor) then
         i_d = 0
         do i = 1, n_descriptors
            do n = lbound(my_descriptor_data%x(i)%ii,1),ubound(my_descriptor_data%x(i)%ii,1)
               i_d = i_d + 1
               if(present(grad_descriptor_index)) grad_descriptor_index(:,i_d) = (/i,my_descriptor_data%x(i)%ii(n)/)
               if(present(grad_descriptor_out)) grad_descriptor_out(:,:,i_d) = my_descriptor_data%x(i)%grad_data(:,:,n)
               if(present(grad_descriptor_pos)) grad_descriptor_pos(:,i_d) = my_descriptor_data%x(i)%pos(:,n)
               if(present(grad_covariance_cutoff)) grad_covariance_cutoff(:,i_d) = my_descriptor_data%x(i)%grad_covariance_cutoff(:,n)
            enddo
         enddo
      endif

      call finalise(my_descriptor_data,error=error)

   endsubroutine descriptor_calc_array

   subroutine bispectrum_SO4_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(bispectrum_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(cplx_2d), dimension(:), allocatable :: U
      type(cplx_3d), dimension(:,:), allocatable :: dU

      complex(dp) :: sub
      complex(dp), dimension(3) :: dsub
      real(dp), dimension(3) :: diff, u_ij
      real(dp) :: r, tmp_cg
      integer :: i, n, n_i, ji, jn, j, m1, m2, j1, j2, m11, m12, m21, m22, i_desc, i_bisp, d, n_descriptors, n_cross, l_n_neighbours
      integer, dimension(3) :: shift
      integer, dimension(116) :: species_map
      logical :: my_do_descriptor, my_do_grad_descriptor

      INIT_ERROR(error)

      call system_timer('bispectrum_SO4_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='bispectrum_SO4_calc args_str')) then
            RAISE_ERROR("bispectrum_SO4_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("bispectrum_SO4_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call cg_initialise(this%j_max, 2)

      call finalise(descriptor_out)

      d = bispectrum_SO4_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif

         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif

      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif

         if(my_do_grad_descriptor) then
            ! dU is not allocated, allocate and zero it
            allocate( dU(0:this%j_max,0:n_neighbours(at,i,max_dist=this%cutoff)) )
            do j = 0, this%j_max
               allocate( dU(j,0)%mm(3,-j:j,-j:j) )
               dU(j,0)%mm = CPLX_ZERO
            enddo

            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            ji = neighbour(at, i, n, jn=jn, distance=r, diff=diff, cosines=u_ij,shift=shift)
            if( r >= this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = ji
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,ji) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif
         enddo

         if(my_do_grad_descriptor) then
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,dU,args_str,error=error)
         else
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,args_str=args_str,error=error)
         endif

         if(my_do_descriptor) then

            i_bisp = 0
            do j1 = 0, this%j_max
               j2 = j1
               !do j2 = 0, this%j_max
                  do j = abs(j1-j2), min(this%j_max,j1+j2)
                     if( mod(j1+j2+j,2) == 1 ) cycle
                     
                     i_bisp = i_bisp + 1

                     !do m1 = -j, j, 2
                     !   do m2 = -j, j, 2
                     !      sub = CPLX_ZERO
                     !      do m11 = max(-j1-m1,-j1), min(j1-m1,j1), 2
                     !         do m21 = max(-j2-m2,-j2), min(j2-m2,j2), 2
                     !            sub = sub + cg_array(j1,m11,j,m1,j1,m11+m1) &
                     !            * cg_array(j2,m21,j,m2,j2,m21+m2) &
                     !            * U(j1)%mm(m11,m11+m1) * U(j2)%mm(m21,m21+m2)
                     !         enddo
                     !      enddo
                     !      descriptor_out%x(i_desc)%data(i_bisp) = descriptor_out%x(i_desc)%data(i_bisp) + sub*conjg(U(j)%mm(-m2,m1))*(-1)**(m2/2)
                     !   enddo
                     !enddo

                     do m1 = -j, j, 2
                        do m2 = -j, j, 2
                           sub = CPLX_ZERO
                           do m11 = max(-j1,m1-j2), min(j1,m1+j2), 2
                              do m12 = max(-j1,m2-j2), min(j1,m2+j2), 2
                                 sub = sub + cg_array(j1,m11,j2,m1-m11,j,m1) &
                                 * cg_array(j1,m12,j2,m2-m12,j,m2) &
                                 * U(j1)%mm(m11,m12) * U(j2)%mm(m1-m11,m2-m12)
                              enddo
                           enddo
                           descriptor_out%x(i_desc)%data(i_bisp) = descriptor_out%x(i_desc)%data(i_bisp) + sub*conjg(U(j)%mm(m1,m2))
                        enddo
                     enddo

                  enddo
               !enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            n_i = 0
            do n = 0, n_neighbours(at,i)
               if( n>0 ) then
                  ji = neighbour(at, i, n, distance=r)
                  if( r >= this%cutoff ) cycle
                  n_i = n_i + 1
               endif
               i_bisp = 0
               do j1 = 0, this%j_max
                  j2 = j1
                  !do j2 = 0, this%j_max
                     do j = abs(j1-j2), min(this%j_max,j1+j2)
                        if( mod(j1+j2+j,2) == 1 ) cycle

                        i_bisp = i_bisp + 1

                        !do m1 = -j, j, 2
                        !   do m2 = -j, j, 2
                        !      sub = CPLX_ZERO
                        !      dsub = CPLX_ZERO

                        !      do m11 = max(-j1-m1,-j1), min(j1-m1,j1), 2
                        !         do m21 = max(-j2-m2,-j2), min(j2-m2,j2), 2
                        !            tmp_cg =  cg_array(j1,m11,j,m1,j1,m11+m1) &
                        !              * cg_array(j2,m21,j,m2,j2,m21+m2)

                        !            sub = sub + tmp_cg &
                        !            * U(j1)%mm(m11,m1+m11) * U(j2)%mm(m21,m2+m21)
                        !            dsub = dsub + tmp_cg &
                        !            * ( dU(j1,n_i)%mm(:,m11,m1+m11) * U(j2)%mm(m21,m2+m21) + &
                        !            U(j1)%mm(m11,m1+m11) * dU(j2,n_i)%mm(:,m21,m2+m21) )
                        !         enddo
                        !      enddo
                        !      descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) = &
                        !      descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) + &
                        !      ( dsub*conjg(U(j)%mm(-m2,m1)) + sub*conjg(dU(j,n_i)%mm(:,-m2,m1)) )*(-1)**(m2/2)
                        !   enddo
                        !enddo
                        do m1 = -j, j, 2
                           do m2 = -j, j, 2
                              sub = CPLX_ZERO
                              dsub = CPLX_ZERO
                              do m11 = max(-j1,m1-j2), min(j1,m1+j2), 2
                                 do m12 = max(-j1,m2-j2), min(j1,m2+j2), 2
                                    
                                    tmp_cg =  cg_array(j1,m11,j2,m1-m11,j,m1) &
                                    * cg_array(j1,m12,j2,m2-m12,j,m2)

                                    sub = sub + tmp_cg &
                                    * U(j1)%mm(m11,m12) * U(j2)%mm(m1-m11,m2-m12)
                                    dsub = dsub + tmp_cg &
                                    * ( dU(j1,n_i)%mm(:,m11,m12) * U(j2)%mm(m1-m11,m2-m12) + &
                                    U(j1)%mm(m11,m12) * dU(j2,n_i)%mm(:,m1-m11,m2-m12) )
                                 enddo
                              enddo
                              descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) = &
                              descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) + &
                              dsub*conjg(U(j)%mm(m1,m2)) + sub*conjg(dU(j,n_i)%mm(:,m1,m2))
                           enddo
                        enddo

                     enddo
                  !enddo
               enddo 
            enddo
         endif

         call finalise(dU)
      enddo ! i

      ! clear U from the memory
      call finalise(U)

      call system_timer('bispectrum_SO4_calc')

   endsubroutine bispectrum_SO4_calc

   subroutine bispectrum_so3_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(bispectrum_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(cplx_1d), dimension(:), allocatable :: SphericalY_ij
      type(cplx_1d), dimension(:,:), allocatable :: fourier_so3

      type(cplx_2d), dimension(:), allocatable :: dSphericalY_ij
      type(cplx_2d), dimension(:,:,:), allocatable :: dfourier_so3

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, a, l, m, l1, l2, m1, i_desc, i_pow, l_n_neighbours, n_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:), allocatable :: Rad_ij
      real(dp), dimension(:,:), allocatable :: dRad_ij

      complex(dp) :: sub, dsub(3)

      integer, dimension(116) :: species_map

      INIT_ERROR(error)

      call system_timer('bispectrum_so3_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_so3_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call cg_initialise(this%l_max)

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='bispectrum_SO3_calc args_str')) then
            RAISE_ERROR("bispectrum_SO3_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("bispectrum_SO3_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("bispectrum_SO3_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = bispectrum_so3_dimensions(this,error)
         
      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif


      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(fourier_so3(0:this%l_max,this%n_max),SphericalY_ij(0:this%l_max),Rad_ij(this%n_max))
      do a = 1, this%n_max
         do l = 0, this%l_max
            allocate(fourier_so3(l,a)%m(-l:l))
            fourier_so3(l,a)%m(:) = CPLX_ZERO
         enddo
      enddo
      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
      enddo

      if(my_do_grad_descriptor) then
         allocate( dRad_ij(3,this%n_max), dSphericalY_ij(0:this%l_max) )
         do l = 0, this%l_max
            allocate(dSphericalY_ij(l)%mm(3,-l:l))
         enddo
      endif

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         do a = 1, this%n_max
            do l = 0, this%l_max
               fourier_so3(l,a)%m(:) = CPLX_ZERO
            enddo
         enddo

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif

         if(my_do_grad_descriptor) then
            allocate( dfourier_so3(0:this%l_max,this%n_max,0:n_neighbours(at,i,max_dist=this%cutoff)) )
            do n = 0, n_neighbours(at,i,max_dist=this%cutoff)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     allocate(dfourier_so3(l,a,n)%mm(3,-l:l))
                     dfourier_so3(l,a,n)%mm(:,:) = CPLX_ZERO
                  enddo
               enddo
            enddo
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a)
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij
            enddo

            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian(l,m,d_ij)
                  if(my_do_grad_descriptor) dSphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian(l,m,d_ij)
               enddo
            enddo

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = -l, l
                     fourier_so3(l,a)%m(m) = fourier_so3(l,a)%m(m) + Rad_ij(a)*SphericalY_ij(l)%m(m)
                     if(my_do_grad_descriptor) then
                        dfourier_so3(l,a,n_i)%mm(:,m) = dfourier_so3(l,a,n_i)%mm(:,m) + &
                        dRad_ij(:,a) * SphericalY_ij(l)%m(m) + Rad_ij(a)*dSphericalY_ij(l)%mm(:,m)
                     endif
                  enddo
               enddo
            enddo

         enddo ! n

         if(my_do_descriptor) then
            i_pow = 0
            do a = 1, this%n_max
               do l1 = 0, this%l_max
                  l2 = l1
                  !do l2 = 0, this%l_max
                     do l = abs(l1-l2), min(this%l_max,l1+l2)
                        if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                        i_pow = i_pow + 1

                        do m = -l, l
                           sub = CPLX_ZERO
                           do m1 = max(-l1,m-l2),min(l1,m+l2)
                              sub = sub + cg_array(l1,m1,l2,m-m1,l,m) * conjg(fourier_so3(l1,a)%m(m1)) * conjg(fourier_so3(l2,a)%m(m-m1))
                           enddo

                           descriptor_out%x(i_desc)%data(i_pow) = descriptor_out%x(i_desc)%data(i_pow) + fourier_so3(l,a)%m(m) * sub
                        enddo

                     enddo
                  !enddo
               enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            do n = 1, n_neighbours(at,i,max_dist=this%cutoff)
               i_pow = 0
               do a = 1, this%n_max
                  do l1 = 0, this%l_max
                     l2 = l1
                     !do l2 = 0, this%l_max
                        do l = abs(l1-l2), min(this%l_max,l1+l2)
                           if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                           i_pow = i_pow + 1

                           do m = -l, l
                              sub = CPLX_ZERO
                              dsub = CPLX_ZERO
                              do m1 = max(-l1,m-l2),min(l1,m+l2)
                                 dsub = dsub + cg_array(l1,m1,l2,m-m1,l,m) * &
                                 ( conjg(dfourier_so3(l1,a,n)%mm(:,m1)) * conjg(fourier_so3(l2,a)%m(m-m1)) + &
                                   conjg(fourier_so3(l1,a)%m(m1)) * conjg(dfourier_so3(l2,a,n)%mm(:,m-m1)) )
                                 sub = sub + cg_array(l1,m1,l2,m-m1,l,m) * conjg(fourier_so3(l1,a)%m(m1)) * conjg(fourier_so3(l2,a)%m(m-m1))
                              enddo

                              descriptor_out%x(i_desc)%grad_data(i_pow,:,n) = descriptor_out%x(i_desc)%grad_data(i_pow,:,n) + &
                              fourier_so3(l,a)%m(m) * dsub + dfourier_so3(l,a,n)%mm(:,m) * sub
                           enddo
                        enddo
                     !enddo
                  enddo
               enddo
               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n)
            enddo
         endif

         if(allocated(dfourier_so3)) then
            do n = lbound(dfourier_so3,3), ubound(dfourier_so3,3)
               do a = lbound(dfourier_so3,2), ubound(dfourier_so3,2)
                  do l = lbound(dfourier_so3,1), ubound(dfourier_so3,1)
                     deallocate(dfourier_so3(l,a,n)%mm)
                  enddo
               enddo
            enddo
            deallocate(dfourier_so3)
         endif

      enddo ! i

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(dRad_ij)) deallocate(dRad_ij)

      if(allocated(fourier_so3)) then
         do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
            do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
               deallocate(fourier_so3(l,a)%m)
            enddo
         enddo
         deallocate(fourier_so3)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(dSphericalY_ij)) then
         do l = lbound(dSphericalY_ij,1), ubound(dSphericalY_ij,1)
            deallocate(dSphericalY_ij(l)%mm)
         enddo
         deallocate(dSphericalY_ij)
      endif

      call system_timer('bispectrum_so3_calc')

   endsubroutine bispectrum_so3_calc

   subroutine behler_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(behler), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, a, b, i_desc, l_n_neighbours, n_i, m_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, Ang, dAng, Rad, dRad_ij, dRad_ik, dRad_jk, f_cut_ij, f_cut_ik, f_cut_jk, df_cut_ij, df_cut_ik, df_cut_jk, g2, dg2
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('behler_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("behler_calc: descriptor object not initialised", error)
      endif

      if( at%cutoff < this%cutoff ) then
         RAISE_ERROR("behler_calc: cutoff of atoms object ("//at%cutoff//") less than cutoff of descriptor ("//this%cutoff//")", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='behler_calc args_str')) then
            RAISE_ERROR("behler_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("behler_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      call finalise(descriptor_out)

      d = behler_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            f_cut_ij = cos_cutoff_function(r_ij,this%cutoff)
            if(my_do_grad_descriptor) df_cut_ij = dcos_cutoff_function(r_ij,this%cutoff)

            do a = 1, this%n_g2
               g2 = exp(-this%g2(a)%eta * (r_ij-this%g2(a)%rs)**2)
               if(my_do_descriptor) descriptor_out%x(i_desc)%data(a) = descriptor_out%x(i_desc)%data(a) + g2 * f_cut_ij
               if(my_do_grad_descriptor) then
                  dg2 = -2.0_dp * this%g2(a)%eta * (r_ij-this%g2(a)%rs) * g2
                  descriptor_out%x(i_desc)%grad_data(a,:,n_i) = ( dg2 * f_cut_ij + g2 * df_cut_ij ) * u_ij
                  descriptor_out%x(i_desc)%grad_data(a,:,0) = descriptor_out%x(i_desc)%grad_data(a,:,0) - descriptor_out%x(i_desc)%grad_data(a,:,n_i)
               endif
            enddo


            m_i = 0
            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik >= this%cutoff ) cycle

               m_i = m_i + 1

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               u_jk = d_jk / r_jk

               f_cut_ik = cos_cutoff_function(r_ik,this%cutoff)
               f_cut_jk = cos_cutoff_function(r_jk,this%cutoff)

               cos_ijk = dot_product(u_ij,u_ik)

               if(my_do_grad_descriptor) then
                  df_cut_ik = dcos_cutoff_function(r_ik,this%cutoff)
                  df_cut_jk = dcos_cutoff_function(r_jk,this%cutoff)

                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik
               endif

               do b = 1, this%n_g3
                  a = b + this%n_g2

                  Ang = (1.0_dp + this%g3(b)%lambda * cos_ijk)**this%g3(b)%zeta
                  Rad = exp( -this%g3(b)%eta * (r_ij**2 + r_ik**2 + r_jk**2) )

                  if(my_do_descriptor) descriptor_out%x(i_desc)%data(a) = descriptor_out%x(i_desc)%data(a) + 0.5_dp * Ang * Rad * f_cut_ij * f_cut_ik * f_cut_jk
                  if(my_do_grad_descriptor) then
                     dAng = this%g3(b)%zeta * (1.0_dp + this%g3(b)%lambda * cos_ijk)**(this%g3(b)%zeta -1.0_dp) * this%g3(b)%lambda
                     dRad_ij = -this%g3(b)%eta * 2.0_dp * r_ij * Rad
                     dRad_ik = -this%g3(b)%eta * 2.0_dp * r_ik * Rad
                     dRad_jk = -this%g3(b)%eta * 2.0_dp * r_jk * Rad

                     descriptor_out%x(i_desc)%grad_data(a,:,n_i) = descriptor_out%x(i_desc)%grad_data(a,:,n_i) + 0.5_dp * &
                     ( ( dAng * dcosijk_ij * Rad + Ang * ( dRad_ij * u_ij - dRad_jk * u_jk ) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_ik * ( df_cut_ij * u_ij * f_cut_jk - f_cut_ij * df_cut_jk * u_jk ) )

                     descriptor_out%x(i_desc)%grad_data(a,:,m_i) = descriptor_out%x(i_desc)%grad_data(a,:,m_i) + 0.5_dp * &
                     ( ( dAng * dcosijk_ik * Rad + Ang * ( dRad_ik * u_ik + dRad_jk * u_jk ) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_ij * ( df_cut_ik * u_ik * f_cut_jk + f_cut_ik * df_cut_jk * u_jk ) ) 

                     descriptor_out%x(i_desc)%grad_data(a,:,0) = descriptor_out%x(i_desc)%grad_data(a,:,0) - 0.5_dp * &
                     ( ( dAng * (dcosijk_ij+dcosijk_ik) * Rad + Ang * (dRad_ij * u_ij + dRad_ik * u_ik) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_jk * ( df_cut_ij * u_ij * f_cut_ik + f_cut_ij * df_cut_ik * u_ik ) )
                  endif


               enddo

            enddo
         enddo

         do b = 1, this%n_g3
            a = b + this%n_g2

            if(my_do_descriptor) descriptor_out%x(i_desc)%data(a) = descriptor_out%x(i_desc)%data(a) * 2.0_dp**(1.0_dp-this%g3(b)%zeta) 
            if(my_do_grad_descriptor) descriptor_out%x(i_desc)%grad_data(a,:,:) = descriptor_out%x(i_desc)%grad_data(a,:,:) * 2.0_dp**(1.0_dp-this%g3(b)%zeta)
         enddo
      enddo

      call system_timer('behler_calc')

   endsubroutine behler_calc

   subroutine distance_2b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical :: needs_resid
      logical, dimension(:), pointer :: atom_mask_pointer
      integer, dimension(:), pointer :: resid_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zi1, Zi2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, n
      integer, dimension(3) :: shift
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij

      INIT_ERROR(error)

      call system_timer('distance_2b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='distance_2b_calc args_str')) then
            RAISE_ERROR("distance_2b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("distance_2b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      needs_resid = this%only_intra .or. this%only_inter
      if (needs_resid) then
         if (.not. assign_pointer(at, trim(this%resid_name), resid_pointer)) then
            RAISE_ERROR("distance_2b_calc did not find "//trim(this%resid_name)//" property (residue id) in the atoms object.", error)
         end if
      else
         resid_pointer => null()
      end if

      d = distance_2b_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(2))
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:1))
            allocate(descriptor_out%x(i)%ii(0:1))
            allocate(descriptor_out%x(i)%pos(3,0:1))
            allocate(descriptor_out%x(i)%has_grad_data(0:1))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,0:1))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)
            if( r_ij >= this%cutoff ) cycle
!if( r_ij < 3.5_dp ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            if (needs_resid) then
               if (this%only_intra .and. resid_pointer(i) /= resid_pointer(j)) cycle
               if (this%only_inter .and. resid_pointer(i) == resid_pointer(j)) cycle
            end if

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%data(1) = r_ij
               descriptor_out%x(i_desc)%ci(1:2) = (/i,j/)
               descriptor_out%x(i_desc)%has_data = .true.

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_ij,this%cutoff,this%cutoff_transition_width) !coordination_function(r_ij, 3.5_dp, 0.5_dp,this%cutoff,0.5_dp) !
            endif
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(0) = i
               descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
               descriptor_out%x(i_desc)%has_grad_data(0) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,0) = -u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = -dcoordination_function(r_ij,this%cutoff,this%cutoff_transition_width)*u_ij !-dcoordination_function(r_ij, 3.5_dp, 0.5_dp,this%cutoff,0.5_dp)*u_ij !

               descriptor_out%x(i_desc)%ii(1) = j
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(1) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0)

            endif
         enddo
      enddo

      call system_timer('distance_2b_calc')

   endsubroutine distance_2b_calc

   subroutine coordination_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(coordination), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, i_n, l_n_neighbours, i_desc, n_descriptors, n_cross
      integer, dimension(3) :: shift
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, df_cut

      INIT_ERROR(error)

      call system_timer('coordination_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='coordination_calc args_str')) then
            RAISE_ERROR("coordination_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("coordination_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      call finalise(descriptor_out)

      d = coordination_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc = i_desc + 1
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%has_data = .false.

            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         i_n = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)

            if( r_ij >= this%cutoff ) cycle
            i_n = i_n + 1

            if(my_do_descriptor) &
               descriptor_out%x(i_desc)%data(1) = descriptor_out%x(i_desc)%data(1) + coordination_function(r_ij,this%cutoff,this%transition_width)

            if(my_do_grad_descriptor) then
               df_cut = dcoordination_function(r_ij,this%cutoff,this%transition_width) * u_ij

               descriptor_out%x(i_desc)%grad_data(1,:,0) = descriptor_out%x(i_desc)%grad_data(1,:,0) - df_cut

               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,i_n) = df_cut
            endif
         enddo
      enddo

      call system_timer('coordination_calc')

   endsubroutine coordination_calc

   subroutine angle_3b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zk1, Zk2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m
      integer, dimension(3) :: shift_ij, shift_ik
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, fc_j, fc_k, dfc_j, dfc_k
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('angle_3b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='angle_3b_calc args_str')) then
            RAISE_ERROR("angle_3b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("angle_3b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      d = angle_3b_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(1))
            descriptor_out%x(i)%has_data = .false.
         endif

         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:2))
            allocate(descriptor_out%x(i)%ii(0:2))
            allocate(descriptor_out%x(i)%pos(3,0:2))
            allocate(descriptor_out%x(i)%has_grad_data(0:2))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,0:2))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, diff = d_ij, shift=shift_ij)

            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            fc_j = coordination_function(r_ij,this%cutoff,0.5_dp)
            dfc_j = dcoordination_function(r_ij,this%cutoff,0.5_dp)

            do m = 1, n_neighbours(at,i)

               if( n == m ) cycle

               k = neighbour(at, i, m, distance = r_ik, cosines = u_ik, diff = d_ik, shift=shift_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               d_jk = d_ij - d_ik
               r_jk = norm(d_jk)
               u_jk = d_jk / r_jk

               fc_k = coordination_function(r_ik,this%cutoff,0.5_dp)
               dfc_k = dcoordination_function(r_ik,this%cutoff,0.5_dp)

               cos_ijk = dot_product(d_ij,d_ik)/(r_ij*r_ik)

               i_desc = i_desc + 1

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(1) = r_ij + r_ik
                  descriptor_out%x(i_desc)%data(2) = (r_ij - r_ik)**2
                  descriptor_out%x(i_desc)%data(3) = r_jk !cos_ijk
                  descriptor_out%x(i_desc)%ci(1) = i
                  descriptor_out%x(i_desc)%has_data = .true.

                  descriptor_out%x(i_desc)%covariance_cutoff = fc_j*fc_k
               endif

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  descriptor_out%x(i_desc)%ii(0) = i
                  descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
                  descriptor_out%x(i_desc)%has_grad_data(0) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,0) = - u_ij - u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,0) = 2.0_dp * (r_ij - r_ik)*(-u_ij + u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,0) = 0.0_dp !-dcosijk_ij - dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = - dfc_j*fc_k*u_ij - dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(1) = j
                  descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift_ij)
                  descriptor_out%x(i_desc)%has_grad_data(1) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij
                  descriptor_out%x(i_desc)%grad_data(2,:,1) = 2.0_dp * (r_ij - r_ik)*u_ij
                  descriptor_out%x(i_desc)%grad_data(3,:,1) = u_jk !dcosijk_ij

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dfc_j*fc_k*u_ij

                  descriptor_out%x(i_desc)%ii(2) = k
                  descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,k) + matmul(at%lattice,shift_ik)
                  descriptor_out%x(i_desc)%has_grad_data(2) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,2) = u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,2) = 2.0_dp * (r_ij - r_ik)*(-u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,2) = -u_jk !dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = dfc_k*fc_j*u_ik
               endif
            enddo
         enddo
      enddo

      call system_timer('angle_3b_calc')

   endsubroutine angle_3b_calc

   subroutine co_angle_3b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(co_angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(descriptor) :: my_coordination
      type(descriptor_data) :: descriptor_coordination

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zk1, Zk2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m, l_n_neighbours_coordination
      integer, dimension(3) :: shift_ij, shift_ik
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, fc_j, fc_k, dfc_j, dfc_k
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('co_angle_3b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='co_angle_3b_calc args_str')) then
            RAISE_ERROR("co_angle_3b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("co_angle_3b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("co_angle_3b_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = co_angle_3b_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         l_n_neighbours_coordination = n_neighbours(at,i,max_dist=this%coordination_cutoff)

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij)
            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, n_neighbours(at,i)

               if( n == m ) cycle

               k = neighbour(at, i, m, distance = r_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               i_desc = i_desc + 1
               if(my_do_descriptor) then
                  allocate(descriptor_out%x(i_desc)%data(d))
                  descriptor_out%x(i_desc)%data = 0.0_dp
                  allocate(descriptor_out%x(i_desc)%ci(1))
                  descriptor_out%x(i_desc)%has_data = .false.
               endif

               if(my_do_grad_descriptor) then

                  allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:2+l_n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%ii(0:2+l_n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%pos(3,0:2+l_n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%has_grad_data(0:2+l_n_neighbours_coordination))
                  descriptor_out%x(i_desc)%grad_data = 0.0_dp
                  descriptor_out%x(i_desc)%ii = 0
                  descriptor_out%x(i_desc)%pos = 0.0_dp
                  descriptor_out%x(i_desc)%has_grad_data = .false.

                  allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:2+l_n_neighbours_coordination))
                  descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
               endif
            enddo
         enddo
      enddo

      call initialise(my_coordination,'coordination cutoff='//this%coordination_cutoff//' coordination_transition_width='//this%coordination_transition_width,error)
      call calc(my_coordination,at,descriptor_coordination,do_descriptor,do_grad_descriptor,args_str,error)
      
      i_desc = 0
      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, diff = d_ij, shift=shift_ij)

            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            fc_j = coordination_function(r_ij,this%cutoff,0.5_dp)
            dfc_j = dcoordination_function(r_ij,this%cutoff,0.5_dp)

            do m = 1, n_neighbours(at,i)
               if( n == m ) cycle

               k = neighbour(at, i, m, distance = r_ik, cosines = u_ik, diff = d_ik, shift=shift_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               d_jk = d_ij - d_ik
               r_jk = norm(d_jk)
               u_jk = d_jk / r_jk

               fc_k = coordination_function(r_ik,this%cutoff,0.5_dp)
               dfc_k = dcoordination_function(r_ik,this%cutoff,0.5_dp)

               cos_ijk = dot_product(d_ij,d_ik)/(r_ij*r_ik)

               i_desc = i_desc + 1

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(1) = r_ij + r_ik
                  descriptor_out%x(i_desc)%data(2) = (r_ij - r_ik)**2
                  descriptor_out%x(i_desc)%data(3) = r_jk !cos_ijk
                  descriptor_out%x(i_desc)%data(4) = descriptor_coordination%x(i)%data(1)
                  descriptor_out%x(i_desc)%ci(1) = i
                  descriptor_out%x(i_desc)%has_data = .true.

                  descriptor_out%x(i_desc)%covariance_cutoff = fc_j*fc_k
               endif

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  descriptor_out%x(i_desc)%ii(0) = i
                  descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
                  descriptor_out%x(i_desc)%has_grad_data(0) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,0) = - u_ij - u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,0) = 2.0_dp * (r_ij - r_ik)*(-u_ij + u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,0) = 0.0_dp !-dcosijk_ij - dcosijk_ik
                  descriptor_out%x(i_desc)%grad_data(4,:,0) = descriptor_coordination%x(i)%grad_data(1,:,0)

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = - dfc_j*fc_k*u_ij - dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(1) = j
                  descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift_ij)
                  descriptor_out%x(i_desc)%has_grad_data(1) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij
                  descriptor_out%x(i_desc)%grad_data(2,:,1) = 2.0_dp * (r_ij - r_ik)*u_ij
                  descriptor_out%x(i_desc)%grad_data(3,:,1) = u_jk !dcosijk_ij

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dfc_j*fc_k*u_ij

                  descriptor_out%x(i_desc)%ii(2) = k
                  descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,k) + matmul(at%lattice,shift_ik)
                  descriptor_out%x(i_desc)%has_grad_data(2) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,2) = u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,2) = 2.0_dp * (r_ij - r_ik)*(-u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,2) = -u_jk !dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(3:) = descriptor_coordination%x(i)%ii(1:)
                  descriptor_out%x(i_desc)%pos(:,3:) = descriptor_coordination%x(i)%pos(:,1:)
                  descriptor_out%x(i_desc)%has_grad_data(3:) = descriptor_coordination%x(i)%has_grad_data(1:)
                  descriptor_out%x(i_desc)%grad_data(4,:,3:) = descriptor_coordination%x(i)%grad_data(1,:,1:)
               endif
            enddo
         enddo
      enddo

      call finalise(my_coordination)
      call finalise(descriptor_coordination)

      call system_timer('co_angle_3b_calc')

   endsubroutine co_angle_3b_calc

   subroutine co_distance_2b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(co_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(descriptor) :: my_coordination
      type(descriptor_data) :: descriptor_coordination

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zi1, Zi2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, n, n_neighbours_coordination_i, n_neighbours_coordination_ij
      integer, dimension(3) :: shift
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij

      INIT_ERROR(error)
      call system_timer('co_distance_2b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='co_distance_2b_calc args_str')) then
            RAISE_ERROR("co_distance_2b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("co_distance_2b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("co_distance_2b_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = co_distance_2b_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         
         if( associated(atom_mask_pointer) ) then
            if( .not. atom_mask_pointer(i) ) cycle
         endif

         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance=r_ij)

            if(r_ij >= this%cutoff) cycle
!if(r_ij <3.5_dp) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc)%data(d))
               descriptor_out%x(i_desc)%data = 0.0_dp
               allocate(descriptor_out%x(i_desc)%ci(2))
               descriptor_out%x(i_desc)%has_data = .false.
            endif

            if(my_do_grad_descriptor) then
               n_neighbours_coordination_ij = n_neighbours(at,i,max_dist=this%coordination_cutoff) + &
               n_neighbours(at,j,max_dist=this%coordination_cutoff) + 2

               allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%ii(0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%pos(3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%has_grad_data(0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_data = 0.0_dp
               descriptor_out%x(i_desc)%ii = 0
               descriptor_out%x(i_desc)%pos = 0.0_dp
               descriptor_out%x(i_desc)%has_grad_data = .false.

               allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
            endif
         enddo
      enddo

      call initialise(my_coordination,'coordination cutoff='//this%coordination_cutoff//' transition_width='//this%coordination_transition_width,error)
      call calc(my_coordination,at,descriptor_coordination,.true.,do_grad_descriptor,args_str,error)
      
      i_desc = 0
      do i = 1, at%N

         if( associated(atom_mask_pointer) ) then
            if( .not. atom_mask_pointer(i) ) cycle
         endif

         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)
            if( r_ij >= this%cutoff ) cycle
!if(r_ij <3.5_dp) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(1:2) = (/i,j/)
               
               descriptor_out%x(i_desc)%has_data = .true.

               descriptor_out%x(i_desc)%data(1) = r_ij
               descriptor_out%x(i_desc)%data(2) = descriptor_coordination%x(i)%data(1) + descriptor_coordination%x(j)%data(1)
               descriptor_out%x(i_desc)%data(3) = (descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))**2

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_ij, this%cutoff,this%transition_width) !coordination_function(r_ij,3.5_dp, 0.5_dp, this%cutoff,this%transition_width)
            endif
            if(my_do_grad_descriptor) then
               n_neighbours_coordination_i = n_neighbours(at,i,max_dist=this%coordination_cutoff)

               descriptor_out%x(i_desc)%ii(0) = i
               descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
               descriptor_out%x(i_desc)%has_grad_data(0) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,0) = -u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = -dcoordination_function(r_ij,this%cutoff,this%transition_width)*u_ij !-dcoordination_function(r_ij,3.5_dp, 0.5_dp, this%cutoff,this%transition_width)*u_ij

               descriptor_out%x(i_desc)%ii(1) = j
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(1) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0)

               descriptor_out%x(i_desc)%ii(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%ii(:)
               descriptor_out%x(i_desc)%pos(:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%pos(:,:)
               descriptor_out%x(i_desc)%has_grad_data(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%has_grad_data(:)
               descriptor_out%x(i_desc)%grad_data(2,:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%grad_data(1,:,:)
               descriptor_out%x(i_desc)%grad_data(3,:,2:n_neighbours_coordination_i+2) = 2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
                  descriptor_coordination%x(i)%grad_data(1,:,:)

               descriptor_out%x(i_desc)%ii(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%ii(:)
               descriptor_out%x(i_desc)%pos(:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%pos(:,:)
               descriptor_out%x(i_desc)%has_grad_data(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%has_grad_data(:)
               descriptor_out%x(i_desc)%grad_data(2,:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%grad_data(1,:,:)
               descriptor_out%x(i_desc)%grad_data(3,:,n_neighbours_coordination_i+3:) = -2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
                  descriptor_coordination%x(j)%grad_data(1,:,:)

            endif
         enddo
      enddo

      call finalise(my_coordination)
      call finalise(descriptor_coordination)

      call system_timer('co_distance_2b_calc')

   endsubroutine co_distance_2b_calc

   subroutine cosnx_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(cosnx), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, a, b, i_desc, i_cosnx, l_n_neighbours, n_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, T_0_cos_ijk, T_1_cos_ijk, T_n_cos_ijk, U_0_cos_ijk, U_1_cos_ijk, U_n_cos_ijk, Ang
      real(dp), dimension(3) :: u_ij, u_ik, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik, dAng_ij, dAng_ik
      real(dp), dimension(:), allocatable :: Rad_ij, Rad_ik, T_cos_ijk, U_cos_ijk
      real(dp), dimension(:,:), allocatable :: dRad_ij, dRad_ik
      integer, dimension(116) :: species_map

      INIT_ERROR(error)

      call system_timer('cosnx_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='cosnx_calc args_str')) then
            RAISE_ERROR("cosnx_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("cosnx_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = cosnx_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(Rad_ij(this%n_max), Rad_ik(this%n_max))
      allocate(T_cos_ijk(0:this%l_max))
      if(my_do_grad_descriptor) then
         allocate(U_cos_ijk(-1:this%l_max))
         allocate(dRad_ij(3,this%n_max), dRad_ik(3,this%n_max))
      endif

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a) * this%w(species_map(at%Z(j)))
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij * this%w(species_map(at%Z(j)))
            enddo

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik >= this%cutoff ) cycle

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik
               endif

               do a = 1, this%n_max
                  Rad_ik(a) = RadialFunction(this%Radial, r_ik, a) * this%w(species_map(at%Z(k)))
                  if(my_do_grad_descriptor) dRad_ik(:,a) = GradRadialFunction(this%Radial, r_ik, a) * u_ik * this%w(species_map(at%Z(k)))
               enddo

               if(this%l_max >= 0) then
                  T_cos_ijk(0) = 1.0_dp
                  T_0_cos_ijk = T_cos_ijk(0)
                  if(my_do_grad_descriptor) then
                     U_cos_ijk(-1) = 0.0_dp
                     U_cos_ijk(0) = 1.0_dp
                     U_0_cos_ijk = U_cos_ijk(0)
                  endif
               endif

               if(this%l_max >= 1) then
                  T_cos_ijk(1) = cos_ijk
                  T_1_cos_ijk = T_cos_ijk(1)
                  if(my_do_grad_descriptor) then
                     U_cos_ijk(1) = 2.0_dp*cos_ijk
                     U_1_cos_ijk = U_cos_ijk(1)
                  endif
               endif

               do b = 2, this%l_max
                  T_n_cos_ijk = 2*cos_ijk*T_1_cos_ijk - T_0_cos_ijk
                  T_0_cos_ijk = T_1_cos_ijk
                  T_1_cos_ijk = T_n_cos_ijk

                  T_cos_ijk(b) = T_n_cos_ijk

                  if(my_do_grad_descriptor) then
                     U_n_cos_ijk = 2*cos_ijk*U_1_cos_ijk - U_0_cos_ijk
                     U_0_cos_ijk = U_1_cos_ijk
                     U_1_cos_ijk = U_n_cos_ijk

                     U_cos_ijk(b) = U_n_cos_ijk
                  endif
               enddo
     
               i_cosnx = 0
               do a = 1, this%n_max
                  do b = 0, this%l_max
                     i_cosnx = i_cosnx + 1

                     Ang = T_cos_ijk(b)

                     if(my_do_descriptor) &
                        descriptor_out%x(i_desc)%data(i_cosnx) = descriptor_out%x(i_desc)%data(i_cosnx) + Rad_ij(a)*Rad_ik(a)*Ang*0.5_dp

                     if(my_do_grad_descriptor) then

                        dAng_ij = b*U_cos_ijk(b-1) * dcosijk_ij
                        dAng_ik = b*U_cos_ijk(b-1) * dcosijk_ik

                        descriptor_out%x(i_desc)%grad_data(i_cosnx,:,0) = descriptor_out%x(i_desc)%grad_data(i_cosnx,:,0) - &
                        ( Rad_ij(a)*Rad_ik(a)*(dAng_ij+dAng_ik) + dRad_ij(:,a)*Rad_ik(a)*Ang + Rad_ij(a)*dRad_ik(:,a)*Ang ) * 0.5_dp

                        descriptor_out%x(i_desc)%grad_data(i_cosnx,:,n_i) = descriptor_out%x(i_desc)%grad_data(i_cosnx,:,n_i) + &
                        (Rad_ij(a)*Rad_ik(a)*dAng_ij + dRad_ij(:,a)*Rad_ik(a)*Ang)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(Rad_ik)) deallocate(Rad_ik)
      if(allocated(T_cos_ijk)) deallocate(T_cos_ijk)
      if(allocated(U_cos_ijk)) deallocate(U_cos_ijk)
      if(allocated(dRad_ij)) deallocate(dRad_ij)
      if(allocated(dRad_ik)) deallocate(dRad_ik)

      call system_timer('cosnx_calc')

   endsubroutine cosnx_calc

   subroutine trihis_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(trihis), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, i_desc
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, Sym_Cor_S, Sym_Cor_A, exp_desc
      real(dp), dimension(3) :: u_ij, u_ik, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik, x, exp_arg, dexp_desc
      real(dp), dimension(3,3) :: dx_j, dx_k

      INIT_ERROR(error)

      call system_timer('trihis_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_calc: descriptor object not initialised", error)
      endif
      RAISE_ERROR("trihis_calc: ab686 noticed that this routine needs updating. Remove this line if you know what you are doing, then proceed.", error)

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='trihis_calc args_str')) then
            RAISE_ERROR("trihis_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("trihis_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("trihis_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = trihis_dimensions(this,error)

      allocate(descriptor_out%x(at%N))
      do i = 1, at%N
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(1))
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%ii(0:n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%pos(3,0:n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%has_grad_data(0:n_neighbours(at,i)))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.
         endif
      enddo

      do i = 1, at%N

         if(my_do_descriptor) then
            descriptor_out%x(i)%ci(1) = i
            descriptor_out%x(i)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i)%ii(0) = i
            descriptor_out%x(i)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i)%has_grad_data(0) = .true.
         endif

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            if(my_do_grad_descriptor) then
               descriptor_out%x(i)%ii(n) = j
               descriptor_out%x(i)%pos(:,n) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i)%has_grad_data(n) = .true.
            endif

            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik >= this%cutoff ) cycle

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               Sym_Cor_S = r_ij + r_ik
               Sym_Cor_A = (r_ij - r_ik)**2

               x = (/Sym_Cor_S, Sym_Cor_A, cos_ijk/)

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  dx_j(:,1) = u_ij
                  dx_j(:,2) = 2.0_dp*(r_ij - r_ik)*u_ij
                  dx_j(:,3) = dcosijk_ij

                  dx_k(:,1) = u_ik
                  dx_k(:,2) = -2.0_dp*(r_ij - r_ik)*u_ik
                  dx_k(:,3) = dcosijk_ik
               endif

               do i_desc = 1, this%n_gauss

                  exp_arg = (x - this%gauss_centre(:,i_desc))/this%gauss_width(:,i_desc)
                  exp_desc = exp(-0.5_dp*sum(exp_arg**2))

                  if(my_do_descriptor) &
                     descriptor_out%x(i)%data(i_desc) = descriptor_out%x(i)%data(i_desc) + exp_desc

                  if(my_do_grad_descriptor) then
                     dexp_desc = -exp_desc * exp_arg / this%gauss_width(:,i_desc)

                     descriptor_out%x(i)%grad_data(i_desc,:,0) = descriptor_out%x(i)%grad_data(i_desc,:,0) - &
                        matmul(dx_j+dx_k,dexp_desc)
                     descriptor_out%x(i)%grad_data(i_desc,:,n) = descriptor_out%x(i)%grad_data(i_desc,:,n) + &
                        2.0_dp*matmul(dx_j,dexp_desc)
                  endif
               enddo
            enddo
         enddo
      enddo

      call system_timer('trihis_calc')

   endsubroutine trihis_calc

   subroutine water_monomer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(water_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, i, iO, iH1, iH2
      integer :: i_desc, mpi_n_procs, mpi_my_proc
      integer, dimension(3) :: shift_1, shift_2
      integer, dimension(:,:), allocatable :: water_monomer_index
      real(dp) :: r1, r2
      real(dp), dimension(3) :: v1, v2, u1, u2

      INIT_ERROR(error)

      call system_timer('water_monomer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='water_monomer_calc args_str')) then
            RAISE_ERROR("water_monomer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("water_monomer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      d = water_monomer_dimensions(this,error)
      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer, error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross, error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(3))
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,3))
            allocate(descriptor_out%x(i)%ii(3))
            allocate(descriptor_out%x(i)%pos(3,3))
            allocate(descriptor_out%x(i)%has_grad_data(3))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,3))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(water_monomer_index(3,count(at%Z==8)))
      call find_water_monomer(at,water_monomer_index,error=error)

      i_desc = 0
      do i = 1, count(at%Z==8)

         iO = water_monomer_index(1,i)
         iH1 = water_monomer_index(2,i)
         iH2 = water_monomer_index(3,i)

         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(iO)) cycle
         endif
         i_desc = i_desc + 1

         v1 = diff_min_image(at,iO,iH1,shift=shift_1)
         v2 = diff_min_image(at,iO,iH2,shift=shift_2)
         r1 = sqrt(dot_product(v1,v1))
         r2 = sqrt(dot_product(v2,v2))
         u1 = v1 / r1
         u2 = v2 / r2

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(:) = water_monomer_index(:,i)
            descriptor_out%x(i_desc)%has_data = .true.
            descriptor_out%x(i_desc)%data(1) = r1+r2 
            descriptor_out%x(i_desc)%data(2) = (r1-r2)**2
            descriptor_out%x(i_desc)%data(3) = dot_product(v1,v2)
         endif

         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(:) = water_monomer_index(:,i)
            descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,iO)
            descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,iH1) + matmul(at%lattice,shift_1)
            descriptor_out%x(i_desc)%pos(:,3) = at%pos(:,iH2) + matmul(at%lattice,shift_2)
            descriptor_out%x(i_desc)%has_grad_data(:) = .true.

            descriptor_out%x(i_desc)%grad_data(1,:,1) = -u1-u2                  ! 1st descriptor wrt rO
            descriptor_out%x(i_desc)%grad_data(1,:,2) =  u1                     ! 1st descriptor wrt rH1
            descriptor_out%x(i_desc)%grad_data(1,:,3) =  u2                     ! 1st descriptor wrt rH2
            descriptor_out%x(i_desc)%grad_data(2,:,1) =  2.0_dp*(r1-r2)*(u2-u1) ! 2nd descriptor wrt rO
            descriptor_out%x(i_desc)%grad_data(2,:,2) =  2.0_dp*(r1-r2)*u1      ! 2nd descriptor wrt rH1
            descriptor_out%x(i_desc)%grad_data(2,:,3) = -2.0_dp*(r1-r2)*u2      ! 2nd descriptor wrt rH2
            descriptor_out%x(i_desc)%grad_data(3,:,1) =  -v1-v2                 ! 3rd descriptor wrt rO
            descriptor_out%x(i_desc)%grad_data(3,:,2) =  v2                     ! 3rd descriptor wrt rH1
            descriptor_out%x(i_desc)%grad_data(3,:,3) =  v1                     ! 3rd descriptor wrt rH2
         endif

      enddo

      deallocate(water_monomer_index)
      call system_timer('water_monomer_calc')

   endsubroutine water_monomer_calc

   subroutine water_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(water_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, n, &
         iAO, iAH1, iAH2, iBO, iBH1, iBH2
      integer :: mpi_n_procs, mpi_my_proc
      integer, dimension(3) :: shift_AO_BO, shift_AO_AH1, shift_AO_AH2, shift_AO_BH1, shift_AO_BH2, &
         shift_BO_AH1, shift_BO_AH2, shift_BO_BH1, shift_BO_BH2, &
         shift_AH1_AH2, shift_AH1_BH1, shift_AH1_BH2, shift_AH2_BH1, shift_AH2_BH2, shift_BH1_BH2
      real(dp), dimension(3) :: diff_AO_BO, diff_AO_AH1, diff_AO_AH2, diff_AO_BH1, diff_AO_BH2, &
         diff_BO_AH1, diff_BO_AH2, diff_BO_BH1, diff_BO_BH2, &
         diff_AH1_AH2, diff_AH1_BH1, diff_AH1_BH2, diff_AH2_BH1, diff_AH2_BH2, diff_BH1_BH2
      integer, dimension(:,:), allocatable :: water_monomer_index
      real(dp) :: r_AO_BO, r_AO_AH1, r_AO_AH2, r_AO_BH1, r_AO_BH2, r_BO_AH1, r_BO_AH2, r_BO_BH1, r_BO_BH2, &
         r_AH1_AH2, r_AH1_BH1, r_AH1_BH2, r_AH2_BH1, r_AH2_BH2, r_BH1_BH2
      integer, dimension(1) :: j_array

      INIT_ERROR(error)

      call system_timer('water_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='water_dimer_calc args_str')) then
            RAISE_ERROR("water_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("water_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      d = water_dimer_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(6))
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,6))
            allocate(descriptor_out%x(i)%ii(6))
            allocate(descriptor_out%x(i)%pos(3,6))
            allocate(descriptor_out%x(i)%has_grad_data(6))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,6))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp

         endif
      enddo

      n_monomers = 0
      do i = 1, at%N
         if(at%Z(i) == 8) n_monomers = n_monomers+1
      enddo

      allocate(water_monomer_index(3,n_monomers))
      call find_water_monomer(at,water_monomer_index,OHH_ordercheck=this%OHH_ordercheck,monomer_cutoff=this%monomer_cutoff,error=error)

      i_desc = 0
      do i = 1, n_monomers
         iAO = water_monomer_index(1,i)
         iAH1 = water_monomer_index(2,i)
         iAH2 = water_monomer_index(3,i)

         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(iAO)) cycle
         endif

         diff_AO_AH1 = diff_min_image(at,iAO,iAH1,shift=shift_AO_AH1)
         diff_AO_AH2 = diff_min_image(at,iAO,iAH2,shift=shift_AO_AH2)
         diff_AH1_AH2 = diff_min_image(at,iAH1,iAH2,shift=shift_AH1_AH2)

         r_AO_AH1 = norm(diff_AO_AH1)
         r_AO_AH2 = norm(diff_AO_AH2)
         r_AH1_AH2 = norm(diff_AH1_AH2)

         do n = 1, n_neighbours(at,iAO)
            iBO = neighbour(at,iAO,n,distance=r_AO_BO, diff=diff_AO_BO, shift=shift_AO_BO )
            if(at%Z(iBO) /= 8) cycle
            if( r_AO_BO >= this%cutoff ) cycle
            i_desc = i_desc + 1
            j_array = find(water_monomer_index(1,:) == iBO)
            j = j_array(1)

            iBH1 = water_monomer_index(2,j)
            iBH2 = water_monomer_index(3,j)

            diff_BO_BH1 = diff_min_image(at,iBO,iBH1,shift=shift_BO_BH1)
            diff_BO_BH2 = diff_min_image(at,iBO,iBH2,shift=shift_BO_BH2)
            diff_BH1_BH2 = diff_min_image(at,iBH1,iBH2,shift=shift_BH1_BH2)

            r_BO_BH1 = norm(diff_BO_BH1)
            r_BO_BH2 = norm(diff_BO_BH2)
            r_BH1_BH2 = norm(diff_BH1_BH2)

            diff_AO_BH1 = diff_AO_BO + diff_BO_BH1
            diff_AO_BH2 = diff_AO_BO + diff_BO_BH2
            shift_AO_BH1 = shift_AO_BO + shift_BO_BH1
            shift_AO_BH2 = shift_AO_BO + shift_BO_BH2

            r_AO_BH1 = norm(diff_AO_BH1)
            r_AO_BH2 = norm(diff_AO_BH2)

            diff_BO_AH1 = -diff_AO_BO + diff_AO_AH1
            diff_BO_AH2 = -diff_AO_BO + diff_AO_AH2

            shift_BO_AH1 = -shift_AO_BO + shift_AO_AH1
            shift_BO_AH2 = -shift_AO_BO + shift_AO_AH2

            r_BO_AH1 = norm(diff_BO_AH1)
            r_BO_AH2 = norm(diff_BO_AH2)

            diff_AH1_BH1 = -diff_AO_AH1 + diff_AO_BO + diff_BO_BH1
            diff_AH1_BH2 = -diff_AO_AH1 + diff_AO_BO + diff_BO_BH2
            diff_AH2_BH1 = -diff_AO_AH2 + diff_AO_BO + diff_BO_BH1
            diff_AH2_BH2 = -diff_AO_AH2 + diff_AO_BO + diff_BO_BH2

            shift_AH1_BH1 = -shift_AO_AH1 + shift_AO_BO + shift_BO_BH1
            shift_AH1_BH2 = -shift_AO_AH1 + shift_AO_BO + shift_BO_BH2
            shift_AH2_BH1 = -shift_AO_AH2 + shift_AO_BO + shift_BO_BH1
            shift_AH2_BH2 = -shift_AO_AH2 + shift_AO_BO + shift_BO_BH2

            r_AH1_BH1 = norm(diff_AH1_BH1)
            r_AH1_BH2 = norm(diff_AH1_BH2)
            r_AH2_BH1 = norm(diff_AH2_BH1)
            r_AH2_BH2 = norm(diff_AH2_BH2)


            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(:) = (/ water_monomer_index(:,i),water_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (/r_AO_BO, &
                  r_AO_AH1, r_AO_AH2, r_AO_BH1, r_AO_BH2, r_BO_AH1, r_BO_AH2, r_BO_BH1, r_BO_BH2, &
                  r_AH1_AH2, r_AH1_BH1, r_AH1_BH2, r_AH2_BH1, r_AH2_BH2, r_BH1_BH2/)

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_AO_BO, &
               this%cutoff,this%cutoff_transition_width)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ water_monomer_index(:,i),water_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,iAO) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,iAH1) + matmul(at%lattice,shift_AO_AH1) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,3) = at%pos(:,iAH2) + matmul(at%lattice,shift_AO_AH2) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,4) = at%pos(:,iBO) + matmul(at%lattice,shift_AO_BO) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,5) = at%pos(:,iBH1) + matmul(at%lattice,shift_AO_BH1) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,6) = at%pos(:,iBH2) + matmul(at%lattice,shift_AO_BH2) ! TODO: Have to figure out how to do this.

               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff_AO_BO / r_AO_BO     ! 1st descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(1,:,4) = -descriptor_out%x(i_desc)%grad_data(1,:,1)        ! 1st descriptor wrt OB

               descriptor_out%x(i_desc)%grad_data(2,:,1) = -diff_AO_AH1 / r_AO_AH1  ! 2nd descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(2,:,2) = -descriptor_out%x(i_desc)%grad_data(2,:,1)        ! 2nd descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff_AO_AH2 / r_AO_AH2  ! 3rd descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)        ! 3rd descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(4,:,1) = -diff_AO_BH1 / r_AO_BH1  ! 4th descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(4,:,5) = -descriptor_out%x(i_desc)%grad_data(4,:,1)        ! 4th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(5,:,1) = -diff_AO_BH2 / r_AO_BH2  ! 5th descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(5,:,6) = -descriptor_out%x(i_desc)%grad_data(5,:,1)        ! 5th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(6,:,4) = -diff_BO_AH1 / r_BO_AH1  ! 6th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -descriptor_out%x(i_desc)%grad_data(6,:,4)        ! 6th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(7,:,4) = -diff_BO_AH2 / r_BO_AH2  ! 7th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(7,:,3) = -descriptor_out%x(i_desc)%grad_data(7,:,4)        ! 7th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(8,:,4) = -diff_BO_BH1 / r_BO_BH1  ! 8th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(8,:,5) = -descriptor_out%x(i_desc)%grad_data(8,:,4)        ! 8th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(9,:,4) = -diff_BO_BH2 / r_BO_BH2  ! 9th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(9,:,6) = -descriptor_out%x(i_desc)%grad_data(9,:,4)        ! 9th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(10,:,2) = -diff_AH1_AH2 / r_AH1_AH2 ! 10th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(10,:,3) = -descriptor_out%x(i_desc)%grad_data(10,:,2)         ! 10th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(11,:,2) = -diff_AH1_BH1 / r_AH1_BH1 ! 11th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(11,:,5) = -descriptor_out%x(i_desc)%grad_data(11,:,2)         ! 11th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(12,:,2) = -diff_AH1_BH2 / r_AH1_BH2 ! 12th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(12,:,6) = -descriptor_out%x(i_desc)%grad_data(12,:,2)         ! 12th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(13,:,3) = -diff_AH2_BH1 / r_AH2_BH1 ! 13th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(13,:,5) = -descriptor_out%x(i_desc)%grad_data(13,:,3)         ! 13th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(14,:,3) = -diff_AH2_BH2 / r_AH2_BH2 ! 14th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(14,:,6) = -descriptor_out%x(i_desc)%grad_data(14,:,3)         ! 14th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(15,:,5) = -diff_BH1_BH2 / r_BH1_BH2 ! 15th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(15,:,6) = -descriptor_out%x(i_desc)%grad_data(15,:,5)         ! 15th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -dcoordination_function(r_AO_BO,&
               this%cutoff,this%cutoff_transition_width) * diff_AO_BO / r_AO_BO
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,4) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1)
            endif
         enddo
      enddo

      deallocate(water_monomer_index)
      call system_timer('water_dimer_calc')

   endsubroutine water_dimer_calc

   subroutine A2_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(A2_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, &
         iA1, iA2, iB1, iB2
      integer, dimension(3) :: shift_A1_A2, shift_A1_B1, shift_A1_B2, shift_A2_B1, shift_A2_B2, shift_B1_B2
      integer, dimension(at%N) :: A2_monomer_index
      real(dp) :: r_A1_A2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2, r_B1_B2

      INIT_ERROR(error)

      call system_timer('A2_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='A2_dimer_calc args_str')) then
            RAISE_ERROR("A2_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("A2_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("A2_dimer_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = A2_dimer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(4))
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,4))
            allocate(descriptor_out%x(i)%ii(4))
            allocate(descriptor_out%x(i)%pos(3,4))
            allocate(descriptor_out%x(i)%has_grad_data(4))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,4))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      n_monomers = count(at%Z == this%atomic_number) / 2

      call find_A2_monomer(at,this%atomic_number, this%monomer_cutoff, A2_monomer_index,error)

      i_desc = 0
      do i = 1, at%N
         iA1 = i
         iA2 = neighbour(at,i,A2_monomer_index(i),distance=r_A1_A2,shift=shift_A1_A2)
         if( iA1 > iA2 ) cycle

         do j = i + 1, at%N
            iB1 = j
            iB2 = neighbour(at,j,A2_monomer_index(j),distance=r_B1_B2,shift=shift_B1_B2)
            if( iB1 > iB2 ) cycle

            r_A1_B1 = distance_min_image(at,iA1,iB1,shift=shift_A1_B1)
            r_A1_B2 = distance_min_image(at,iA1,iB2,shift=shift_A1_B2)

            r_A2_B1 = distance_min_image(at,iA2,iB1,shift=shift_A2_B1)
            r_A2_B2 = distance_min_image(at,iA2,iB2,shift=shift_A2_B2)
            
            if( any( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) >= this%cutoff) ) cycle
            i_desc = i_desc + 1

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(:) = (/ iA1, iA2, iB1, iB2 /)
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (/ r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2/)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ iA1, iA2, iB1, iB2 /)
               descriptor_out%x(i_desc)%pos(:,:) = 0.0_dp ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff(at,iA1,iA2,shift=shift_A1_A2) / r_A1_A2      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(1,:,2) = -descriptor_out%x(i_desc)%grad_data(1,:,1)         ! 1st descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(2,:,3) = -diff(at,iB1,iB2,shift=shift_B1_B2) / r_B1_B2      ! 2nd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(2,:,4) = -descriptor_out%x(i_desc)%grad_data(2,:,3)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff(at,iA1,iB1,shift=shift_A1_B1) / r_A1_B1      ! 3rd descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)         ! 3rd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(4,:,1) = -diff(at,iA1,iB2,shift=shift_A1_B2) / r_A1_B2      ! 4th descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(4,:,4) = -descriptor_out%x(i_desc)%grad_data(4,:,1)         ! 4th descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(5,:,2) = -diff(at,iA2,iB1,shift=shift_A2_B1) / r_A2_B1      ! 5th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(5,:,3) = -descriptor_out%x(i_desc)%grad_data(5,:,2)         ! 5th descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -diff(at,iA2,iB2,shift=shift_A2_B2) / r_A2_B2      ! 6th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(6,:,4) = -descriptor_out%x(i_desc)%grad_data(6,:,2)         ! 6th descriptor wrt B2

            endif
         enddo
      enddo

      call system_timer('A2_dimer_calc')

   endsubroutine A2_dimer_calc

   subroutine AB_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(AB_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, &
         iA1, iA2, iB1, iB2
      integer, dimension(3) :: shift_A1_A2, shift_A1_B1, shift_A1_B2, shift_A2_B1, shift_A2_B2, shift_B1_B2
      integer, dimension(:,:), allocatable :: AB_monomer_index
      real(dp) :: r_A1_A2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2, r_B1_B2

      INIT_ERROR(error)

      call system_timer('AB_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='AB_dimer_calc args_str')) then
            RAISE_ERROR("AB_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("AB_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("AB_dimer_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = AB_dimer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(4))
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,4))
            allocate(descriptor_out%x(i)%ii(4))
            allocate(descriptor_out%x(i)%pos(3,4))
            allocate(descriptor_out%x(i)%has_grad_data(4))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,4))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      if( count(at%Z == this%atomic_number1) == count(at%Z == this%atomic_number2) ) then
         n_monomers = count(at%Z == this%atomic_number1)
      else
         RAISE_ERROR("AB_dimer_calc: number of monomer atoms 1 ("//count(at%Z == this%atomic_number1)//") not equal to number of monomer atoms 2 ("//count(at%Z == this%atomic_number1)//")",error)
      endif

      allocate(AB_monomer_index(2,n_monomers))
      call find_AB_monomer(at,(/this%atomic_number1,this%atomic_number2/), this%monomer_cutoff, AB_monomer_index,error)

      i_desc = 0
      do i = 1, n_monomers
         iA1 = AB_monomer_index(1,i)
         iB1 = AB_monomer_index(2,i)
         do j = i + 1, n_monomers
            iA2 = AB_monomer_index(1,j)
            iB2 = AB_monomer_index(2,j)


            r_A1_B1 = distance_min_image(at,iA1,iB1,shift=shift_A1_B1)
            r_A2_B2 = distance_min_image(at,iA2,iB2,shift=shift_A2_B2)

            r_A1_A2 = distance_min_image(at,iA1,iA2,shift=shift_A1_A2)
            r_B1_B2 = distance_min_image(at,iB1,iB2,shift=shift_B1_B2)

            r_A1_B2 = distance_min_image(at,iA1,iB2,shift=shift_A1_B2)
            r_A2_B1 = distance_min_image(at,iA2,iB1,shift=shift_A2_B1)
            
            if( any( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) >= this%cutoff) ) cycle
            i_desc = i_desc + 1

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(:) = (/ AB_monomer_index(:,i),AB_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (/ r_A1_B1, r_A2_B2, r_A1_A2, r_B1_B2, r_A1_B2, r_A2_B1 /)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ AB_monomer_index(:,i),AB_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%pos(:,:) = 0.0_dp ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff(at,iA1,iB1,shift=shift_A1_B1) / r_A1_B1      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(1,:,2) = -descriptor_out%x(i_desc)%grad_data(1,:,1)         ! 1st descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(2,:,3) = -diff(at,iA2,iB2,shift=shift_A2_B2) / r_A2_B2      ! 2nd descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(2,:,4) = -descriptor_out%x(i_desc)%grad_data(2,:,3)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff(at,iA1,iA2,shift=shift_A1_A2) / r_A1_A2      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)         ! 1st descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(4,:,2) = -diff(at,iB1,iB2,shift=shift_B1_B2) / r_B1_B2      ! 2nd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(4,:,4) = -descriptor_out%x(i_desc)%grad_data(4,:,2)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(5,:,1) = -diff(at,iA1,iB2,shift=shift_A1_B2) / r_A1_B2      ! 4th descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(5,:,4) = -descriptor_out%x(i_desc)%grad_data(5,:,1)         ! 4th descriptor wrt B2
               descriptor_out%x(i_desc)%grad_data(6,:,3) = -diff(at,iA2,iB1,shift=shift_A2_B1) / r_A2_B1      ! 5th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -descriptor_out%x(i_desc)%grad_data(6,:,3)         ! 5th descriptor wrt B1

            endif
         enddo
      enddo

      deallocate(AB_monomer_index)
      call system_timer('AB_dimer_calc')

   endsubroutine AB_dimer_calc

   subroutine bond_real_space_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(bond_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(atoms) :: at_copy
      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: n_descriptors, n_cross, i_desc, i, j, n, k, m, m_index, l, ij_neighbours
      integer, dimension(3) :: shift_j, shift_k
      real(dp) :: r_ij, r_ijk
      real(dp) :: atom_i(3), atom_j(3), atom_k(3), bond(3), bond_len
      real(dp) :: atom_i_cross_atom_j(3), atom_i_normsq_min_atom_j_normsq
      real(dp), allocatable :: r(:,:), z(:), c(:)
      real(dp) :: self_overlap
      real(dp), allocatable :: dr(:,:,:,:), dz(:,:,:), dc(:,:,:)
      real(dp), allocatable :: dself_overlap(:,:)
      integer, allocatable :: ii(:)
      real(dp), allocatable :: pos(:,:)
      real(dp) :: r_m_cross_r_l(3)

      INIT_ERROR(error)

      call system_timer('bond_real_space_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='bond_real_space_calc args_str')) then
            RAISE_ERROR("bond_real_space_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("bond_real_space_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("bond_real_space_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0

      do i = 1, at%N
         do n = 1, n_neighbours(at, i)
            j = neighbour(at, i, n, shift=shift_j, distance=r_ij, max_dist=this%bond_cutoff)

            if(j == 0) cycle

            i_desc = i_desc + 1

            atom_i = at%pos(:,i)
            atom_j = at%pos(:,j) + matmul(at%lattice, shift_j)

            at_copy = at
            call add_atoms(at_copy, 0.5_dp * (atom_i + atom_j), 1)
            call calc_connect(at_copy)

            ij_neighbours = 0

            do m = 1, n_neighbours(at_copy, at%N + 1)
               k = neighbour(at_copy, at%N + 1, m, max_dist=this%cutoff)

               if(k == 0) cycle

               if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

               ij_neighbours = ij_neighbours + 1
            enddo

            if(ij_neighbours > this%max_neighbours) then
               RAISE_ERROR("bond_real_space_calc: number of neighbours exceeds max_neighbours", error)
            endif

            if(my_do_descriptor .or. my_do_grad_descriptor) then
               allocate(r(3,ij_neighbours), z(ij_neighbours), c(ij_neighbours))
               allocate(ii(ij_neighbours), pos(3,ij_neighbours))

               r = 0.0_dp
               z = 0.0_dp
               c = 0.0_dp
               self_overlap = 0.0_dp
               bond = atom_i - atom_j
               bond_len = norm(bond)
               atom_i_cross_atom_j = atom_i .cross. atom_j
               atom_i_normsq_min_atom_j_normsq = normsq(atom_i) - normsq(atom_j)
               ii = 0
               pos = 0.0_dp

               if(my_do_grad_descriptor) then
                  allocate(dr(3,ij_neighbours,3,ij_neighbours), dz(ij_neighbours,3,ij_neighbours), dc(ij_neighbours,3,ij_neighbours), dself_overlap(3,ij_neighbours))

                  dr = 0.0_dp
                  dz = 0.0_dp
                  dc = 0.0_dp
                  dself_overlap = 0.0_dp
               endif

               m_index = 2

               do m = 1, n_neighbours(at_copy, at%N + 1)
                  k = neighbour(at_copy, at%N + 1, m, shift=shift_k, distance=r_ijk, max_dist=this%cutoff)

                  if(k == 0) cycle

                  if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

                  atom_k = at_copy%pos(:,k) + matmul(at_copy%lattice, shift_k)

                  if(atom_k .feq. atom_i) then
                     ! r remains zero
                     z(1) = 0.5_dp * bond_len
                     c(1) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     ii(1) = k
                     pos(:,1) = atom_k

                     if(my_do_grad_descriptor) then
                        ! dr remains zero
                        dz(1,:,1) = 0.5_dp * bond / bond_len
                        dz(1,:,2) = - dz(1,:,1)
                        dc(1,:,1) = 0.25_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * bond / r_ijk
                        dc(1,:,2) = - dc(1,:,1)
                     endif
                  elseif(atom_k .feq. atom_j) then
                     ! r remain zero
                     z(2) = -0.5_dp * bond_len
                     c(2) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     ii(2) = k
                     pos(:,2) = atom_k

                     if(my_do_grad_descriptor) then
                        ! dr remains zero
                        dz(2,:,1) = -0.5_dp * bond / bond_len
                        dz(2,:,2) = - dz(2,:,1)
                        dc(2,:,1) = -0.25_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * bond / r_ijk
                        dc(2,:,2) = - dc(2,:,1)
                     endif
                  else
                     m_index = m_index + 1

                     r(:,m_index) = ((atom_k .cross. bond) + atom_i_cross_atom_j) / bond_len
                     z(m_index) = ((atom_k .dot. bond) - 0.5_dp * atom_i_normsq_min_atom_j_normsq) / bond_len
                     c(m_index) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     ii(m_index) = k
                     pos(:,m_index) = atom_k

                     if(my_do_grad_descriptor) then
                        dr(:,m_index,1,1) = ((/ 0.0_dp, atom_k(3) - atom_j(3), atom_j(2) - atom_k(2) /) / bond_len) - (r(:,m_index) * bond(1) / bond_len**2)
                        dr(:,m_index,2,1) = ((/ atom_j(3) - atom_k(3), 0.0_dp, atom_k(1) - atom_j(1) /) / bond_len) - (r(:,m_index) * bond(2) / bond_len**2)
                        dr(:,m_index,3,1) = ((/ atom_k(2) - atom_j(2), atom_j(1) - atom_k(1), 0.0_dp /) / bond_len) - (r(:,m_index) * bond(3) / bond_len**2)
                        dz(m_index,:,1) = ((atom_k - atom_i) / bond_len) - (z(m_index) * bond / bond_len**2)
                        dc(m_index,:,1) = -0.5_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * (atom_k - at_copy%pos(:,at%N + 1)) / r_ijk

                        dr(:,m_index,1,2) = - dr(:,m_index,1,1) + ((/ 0.0_dp, bond(3), - bond(2) /) / bond_len)
                        dr(:,m_index,2,2) = - dr(:,m_index,2,1) + ((/ - bond(3), 0.0_dp, bond(1) /) / bond_len)
                        dr(:,m_index,3,2) = - dr(:,m_index,3,1) + ((/ bond(2), - bond(1), 0.0_dp /) / bond_len)
                        dz(m_index,:,2) = - dz(m_index,:,1) - (bond / bond_len)
                        dc(m_index,:,2) = dc(m_index,:,1)

                        dr(:,m_index,1,m_index) = (/ 0.0_dp, - bond(3), bond(2) /) / bond_len
                        dr(:,m_index,2,m_index) = (/ bond(3), 0.0_dp, - bond(1) /) / bond_len
                        dr(:,m_index,3,m_index) = (/ - bond(2), bond(1), 0.0_dp /) / bond_len
                        dz(m_index,:,m_index) = bond / bond_len
                        dc(m_index,:,m_index) = -2.0_dp * dc(m_index,:,1)
                     endif
                  endif
               enddo
            endif

            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc)%data(2 + (1 + 2 * this%max_neighbours) * this%max_neighbours))
               allocate(descriptor_out%x(i_desc)%ci(2))

               descriptor_out%x(i_desc)%data = 0.0_dp

               do m = 1, ij_neighbours
                  self_overlap = self_overlap + c(m)**2

                  if(m == 1) then
                     descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m)) = z(m)
                  elseif(m == 2) then
                     descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m)) = z(m)

                     self_overlap = self_overlap + 2.0_dp * c(m) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 )
                  else
                     do l = 3, ij_neighbours
                        if(l == m) then
                           descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1) = normsq(r(:,m))
                           descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l)) = z(m)
                        else
                           descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1) = r(:,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l)) = ((r(:,m) .cross. r(:,l)) .dot. bond) / bond_len
                        endif

                        if(l < m) then
                           self_overlap = self_overlap + 2.0_dp * c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                        endif
                     enddo
                  endif
               enddo

               descriptor_out%x(i_desc)%data(1) = real(ij_neighbours, dp)

               descriptor_out%x(i_desc)%data(2) = self_overlap

               descriptor_out%x(i_desc)%data(3:ij_neighbours + 2) = c

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_ij, this%bond_cutoff, this%bond_transition_width)

               descriptor_out%x(i_desc)%ci(:) = (/ i, j /)
               descriptor_out%x(i_desc)%has_data = .true.
            endif

            if(my_do_grad_descriptor) then
               allocate(descriptor_out%x(i_desc)%grad_data(2 + (1 + 2 * ij_neighbours) * ij_neighbours,3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%ii(ij_neighbours))
               allocate(descriptor_out%x(i_desc)%pos(3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%has_grad_data(ij_neighbours))

               descriptor_out%x(i_desc)%grad_data = 0.0_dp

               do m = 1, ij_neighbours
                  dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * c(m) * dc(m,:,1)
                  dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * c(m) * dc(m,:,2)

                  if(m == 1) then
                     descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,1) = dz(m,:,1)
                     descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,2) = dz(m,:,2)
                  elseif(m == 2) then
                     descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,1) = dz(m,:,1)
                     descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,2) = dz(m,:,2)

                     dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * dc(m,:,1) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             + 2.0_dp * c(m) * dc(m - 1,:,1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             + c(m) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             * (pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2
                     dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * dc(m,:,2) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             + 2.0_dp * c(m) * dc(m - 1,:,2) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             + c(m) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             * (pos(:,m - 1) - pos(:,m)) / this%atom_sigma**2
                  else
                     dself_overlap(:,m) = dself_overlap(:,m) + 2.0_dp * c(m) * dc(m,:,m)

                     do l = 3, ij_neighbours
                        if(l == m) then
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,1) = dz(m,:,1)

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,2) = dz(m,:,2)

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,m))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,m))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,m))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,m) = dz(m,:,m)
                        else
                           r_m_cross_r_l = r(:,m) .cross. r(:,l)

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,1) = (dr(:,m,1,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,1,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,1) = (dr(:,m,2,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,2,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,1) = (dr(:,m,3,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,3,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,1) = ((((dr(:,m,1,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,1,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 1.0_dp, 0.0_dp, 0.0_dp /) - (bond * bond(1) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,1) = ((((dr(:,m,2,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,2,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 1.0_dp, 0.0_dp /) - (bond * bond(2) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,1) = ((((dr(:,m,3,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,3,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 0.0_dp, 1.0_dp /) - (bond * bond(3) / bond_len**2)))) / bond_len

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,2) = (dr(:,m,1,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,1,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,2) = (dr(:,m,2,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,2,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,2) = (dr(:,m,3,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,3,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,2) = ((((dr(:,m,1,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,1,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ -1.0_dp, 0.0_dp, 0.0_dp /) + (bond * bond(1) / bond_len**2)))) / bond_len 
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,2) = ((((dr(:,m,2,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,2,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, -1.0_dp, 0.0_dp /) + (bond * bond(2) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,2) = ((((dr(:,m,3,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,3,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 0.0_dp, -1.0_dp /) + (bond * bond(3) / bond_len**2)))) / bond_len

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,m) = dr(:,m,1,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,m) = dr(:,m,2,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,m) = dr(:,m,3,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,m) = ((dr(:,m,1,m) .cross. r(:,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,m) = ((dr(:,m,2,m) .cross. r(:,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,m) = ((dr(:,m,3,m) .cross. r(:,l)) .dot. bond) / bond_len

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,l) = r(:,m) .dot. dr(:,l,1,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,l) = r(:,m) .dot. dr(:,l,2,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,l) = r(:,m) .dot. dr(:,l,3,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,l) = ((r(:,m) .cross. dr(:,l,1,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,l) = ((r(:,m) .cross. dr(:,l,2,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,l) = ((r(:,m) .cross. dr(:,l,3,l)) .dot. bond) / bond_len
                        endif

                        if(l < m) then
                           dself_overlap(:,m) = dself_overlap(:,m) + 2.0_dp * dc(m,:,m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                   + c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                   * (pos(:,l) - pos(:,m)) / this%atom_sigma**2

                           if(l == 1) then
                              dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * dc(m,:,1) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      * (pos(:,m) - pos(:,l)) / this%atom_sigma**2
                              dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * dc(m,:,2) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,2) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                           elseif(l == 2) then
                              dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * dc(m,:,1) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                              dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * dc(m,:,2) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,2) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      * (pos(:,m) - pos(:,l)) / this%atom_sigma**2
                           else
                              dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * dc(m,:,1) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                              dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * dc(m,:,2) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,2) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                              dself_overlap(:,l) = dself_overlap(:,l) + 2.0_dp * c(m) * dc(l,:,l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      * (pos(:,m) - pos(:,l)) / this%atom_sigma**2
                           endif
                        endif
                     enddo
                  endif
               enddo

               !descriptor_out%x(i_desc)%grad_data(1,:,:) = 0.0_dp

               descriptor_out%x(i_desc)%grad_data(2,:,:) = dself_overlap

               descriptor_out%x(i_desc)%grad_data(3:ij_neighbours + 2,:,:) = dc

               descriptor_out%x(i_desc)%ii = ii
               descriptor_out%x(i_desc)%pos = pos

               descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp

               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dcoordination_function(r_ij, this%bond_cutoff, this%bond_transition_width) * bond / r_ij
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = - descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1)

               descriptor_out%x(i_desc)%has_grad_data = .true.
            endif

            if(my_do_descriptor .or. my_do_grad_descriptor) then
               deallocate(r, z, c)
               deallocate(ii, pos)

               if(my_do_grad_descriptor) then
                  deallocate(dr, dz, dc, dself_overlap)
               endif
            endif

            call finalise(at_copy)
         enddo
      enddo

      call system_timer('bond_real_space_calc')

   endsubroutine bond_real_space_calc

   subroutine atom_real_space_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(atom_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, grad_d, n_descriptors, n_cross, descriptor_mould_size, i_desc, i_data, i, j, k, n, l, m, l_n_neighbours, i_n

      real(dp) :: r
      real(dp), dimension(3) :: diff
      real(dp), dimension(1) :: descriptor_mould
      integer, dimension(3) :: shift

      complex(dp), dimension(:), allocatable :: spherical_harmonics
      complex(dp), dimension(:,:), allocatable :: grad_spherical_harmonics

      INIT_ERROR(error)

      call system_timer('atom_real_space_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='atom_real_space_calc args_str')) then
            RAISE_ERROR("atom_real_space_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("atom_real_space_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("atom_real_space_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1

         l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)
         d = ( 2 * (this%l_max+1)**2 + 2 ) * l_n_neighbours

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            grad_d = 2 * (this%l_max+1)**2 + 2

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,1:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(1:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,1:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(1:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,1:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(spherical_harmonics(-this%l_max:this%l_max))
      if( my_do_grad_descriptor ) allocate(grad_spherical_harmonics(3,-this%l_max:this%l_max))

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1
         i_data = 0
         i_n = 0

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif

         if(my_do_grad_descriptor) then
            !descriptor_out%x(i_desc)%ii(0) = i
            !descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            !descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         do n = 1, n_neighbours(at,i)

            j = neighbour(at,i,n,distance = r, diff = diff, shift=shift)
            if(r >= this%cutoff) cycle
            i_n = i_n + 1

            i_data = i_data + 1
            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%data(i_data) = r
            endif
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.
               descriptor_out%x(i_desc)%grad_data(i_data,:,i_n) = diff / r
            endif

            i_data = i_data + 1
            if(my_do_descriptor) descriptor_out%x(i_desc)%data(i_data) = real(i_n,dp)
            if(my_do_grad_descriptor) descriptor_out%x(i_desc)%grad_data(i_data,:,i_n) = real(i_n,dp)

            do l = 0, this%l_max
               descriptor_mould_size = size(transfer(spherical_harmonics(-l:l),descriptor_mould))
               
               do m = -l, l
                  if(my_do_descriptor) spherical_harmonics(m) = SphericalYCartesian(l,m,diff)
                  if(my_do_grad_descriptor) grad_spherical_harmonics(:,m) = GradSphericalYCartesian(l,m,diff)
               enddo

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(i_data+1:i_data+descriptor_mould_size) = transfer(spherical_harmonics(-l:l),descriptor_mould)
               endif

               if(my_do_grad_descriptor) then
                  do k = 1, 3
                     descriptor_out%x(i_desc)%grad_data(i_data+1:i_data+descriptor_mould_size,k,i_n) = &
                     transfer(grad_spherical_harmonics(k,-l:l),descriptor_mould)
                  enddo
               endif

               i_data = i_data + descriptor_mould_size

            enddo
         enddo
      enddo

      if(allocated(spherical_harmonics)) deallocate(spherical_harmonics)
      if(allocated(grad_spherical_harmonics)) deallocate(grad_spherical_harmonics)

      call system_timer('atom_real_space_calc')

   endsubroutine atom_real_space_calc

   subroutine power_so3_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(power_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(cplx_1d), dimension(:), allocatable :: SphericalY_ij
      type(cplx_1d), dimension(:,:), allocatable :: fourier_so3

      type(cplx_2d), dimension(:), allocatable :: dSphericalY_ij
      type(cplx_2d), dimension(:,:,:), allocatable :: dfourier_so3

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, a, l, m, i_desc, i_pow, l_n_neighbours, n_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:), allocatable :: Rad_ij
      real(dp), dimension(:,:), allocatable :: dRad_ij
      integer, dimension(116) :: species_map

      INIT_ERROR(error)

      call system_timer('power_so3_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='power_so3_calc args_str')) then
            RAISE_ERROR("power_so3_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("power_so3_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = power_so3_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(fourier_so3(0:this%l_max,this%n_max),SphericalY_ij(0:this%l_max),Rad_ij(this%n_max))
      do a = 1, this%n_max
         do l = 0, this%l_max
            allocate(fourier_so3(l,a)%m(-l:l))
            fourier_so3(l,a)%m(:) = CPLX_ZERO
         enddo
      enddo
      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
      enddo

      if(my_do_grad_descriptor) then
         allocate( dRad_ij(3,this%n_max), dSphericalY_ij(0:this%l_max) )
         do l = 0, this%l_max
            allocate(dSphericalY_ij(l)%mm(3,-l:l))
         enddo
      endif

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         do a = 1, this%n_max
            do l = 0, this%l_max
               fourier_so3(l,a)%m(:) = CPLX_ZERO
            enddo
         enddo

         if(my_do_grad_descriptor) then
            allocate( dfourier_so3(0:this%l_max,this%n_max,0:n_neighbours(at,i,max_dist=this%cutoff)) )
            do n = 0, n_neighbours(at,i,max_dist=this%cutoff)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     allocate(dfourier_so3(l,a,n)%mm(3,-l:l))
                     dfourier_so3(l,a,n)%mm(:,:) = CPLX_ZERO
                  enddo
               enddo
            enddo
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a)
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij
            enddo

            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian(l,m,d_ij)
                  if(my_do_grad_descriptor) dSphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian(l,m,d_ij)
               enddo
            enddo

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = -l, l
                     fourier_so3(l,a)%m(m) = fourier_so3(l,a)%m(m) + Rad_ij(a)*SphericalY_ij(l)%m(m)
                     if(my_do_grad_descriptor) then
                        dfourier_so3(l,a,n_i)%mm(:,m) = dfourier_so3(l,a,n_i)%mm(:,m) + &
                        dRad_ij(:,a) * SphericalY_ij(l)%m(m) + Rad_ij(a)*dSphericalY_ij(l)%mm(:,m)
                     endif
                  enddo
               enddo
            enddo

         enddo ! n

         if(my_do_descriptor) then
            i_pow = 0
            do a = 1, this%n_max
               do l = 0, this%l_max
                  i_pow = i_pow + 1

                  descriptor_out%x(i_desc)%data(i_pow) = dot_product(fourier_so3(l,a)%m,fourier_so3(l,a)%m)
               enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            do n = 1, n_neighbours(at,i,max_dist=this%cutoff)
               i_pow = 0
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     i_pow = i_pow + 1

                     descriptor_out%x(i_desc)%grad_data(i_pow,:,n) = 2.0_dp * matmul(conjg(dfourier_so3(l,a,n)%mm(:,:)),fourier_so3(l,a)%m(:))
                  enddo
               enddo
               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n)
            enddo
         endif

         if(allocated(dfourier_so3)) then
            do n = lbound(dfourier_so3,3), ubound(dfourier_so3,3)
               do a = lbound(dfourier_so3,2), ubound(dfourier_so3,2)
                  do l = lbound(dfourier_so3,1), ubound(dfourier_so3,1)
                     deallocate(dfourier_so3(l,a,n)%mm)
                  enddo
               enddo
            enddo
            deallocate(dfourier_so3)
         endif

      enddo ! i

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(dRad_ij)) deallocate(dRad_ij)

      if(allocated(fourier_so3)) then
         do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
            do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
               deallocate(fourier_so3(l,a)%m)
            enddo
         enddo
         deallocate(fourier_so3)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(dSphericalY_ij)) then
         do l = lbound(dSphericalY_ij,1), ubound(dSphericalY_ij,1)
            deallocate(dSphericalY_ij(l)%mm)
         enddo
         deallocate(dSphericalY_ij)
      endif

      call system_timer('power_so3_calc')

   endsubroutine power_so3_calc

   subroutine power_SO4_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(power_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(cplx_2d), dimension(:), allocatable :: U
      type(cplx_3d), dimension(:,:), allocatable :: dU

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      real(dp), dimension(3) :: diff, u_ij
      real(dp) :: r
      integer :: i, n, n_i, ji, jn, k, j, i_desc, i_bisp, d, n_descriptors, n_cross, l_n_neighbours
      integer, dimension(3) :: shift
      integer, dimension(116) :: species_map
      logical :: my_do_descriptor, my_do_grad_descriptor

      INIT_ERROR(error)

      call system_timer('power_SO4_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='power_SO4_calc args_str')) then
            RAISE_ERROR("power_SO4_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("power_SO4_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("power_SO4_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = power_SO4_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle

         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif

         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif

      enddo

      i_desc = 0
      do i = 1, at%N

         if( associated(atom_mask_pointer) ) then
            if( .not. atom_mask_pointer(i) ) cycle
         endif

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            ji = neighbour(at, i, n, jn=jn, distance=r, diff=diff, cosines=u_ij,shift=shift)
            if( r >= this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = ji
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,ji) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif
         enddo

         if(my_do_grad_descriptor) then
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,dU,args_str,error=error)
         else
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,args_str=args_str,error=error)
         endif

         if(my_do_descriptor) then

            i_bisp = 0
            do j = 0, this%j_max
               i_bisp = i_bisp + 1
               descriptor_out%x(i_desc)%data(i_bisp) =  sum( conjg(U(j)%mm)*U(j)%mm )
            enddo
         endif

         if(my_do_grad_descriptor) then
            n_i = 0
            do n = 1, n_neighbours(at,i)
               ji = neighbour(at, i, n, distance=r)
               if( r >= this%cutoff ) cycle
               n_i = n_i + 1
               i_bisp = 0
               do j = 0, this%j_max
                  i_bisp = i_bisp + 1
                  do k = 1, 3
                     descriptor_out%x(i_desc)%grad_data(i_bisp,k,n_i) = 2.0_dp * sum( conjg(U(j)%mm)*dU(j,n_i)%mm(k,:,:) )
                  enddo
               enddo 
            enddo
            descriptor_out%x(i_desc)%grad_data(:,:,0) = -sum(descriptor_out%x(i_desc)%grad_data(:,:,:), dim=3)
         endif

         call finalise(dU)
      enddo ! i

      ! clear U from the memory
      call finalise(U)

      call system_timer('power_SO4_calc')

   endsubroutine power_SO4_calc

   subroutine soap_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type real_2d_array
         type(real_2d), dimension(:,:,:), allocatable :: x
      endtype real_2d_array

      type(soap), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(cplx_1d), dimension(:), allocatable, save :: SphericalY_ij
      type(cplx_2d), dimension(:), allocatable, save :: grad_SphericalY_ij

      !SPEED type(cplx_1d), dimension(:,:,:), allocatable :: fourier_so3
      !SPEED type(cplx_2d), dimension(:,:,:), allocatable :: grad_fourier_so3
      type(real_1d), dimension(:,:,:), allocatable, save :: fourier_so3_r, fourier_so3_i, global_fourier_so3_r, global_fourier_so3_i
      type(real_2d), dimension(:,:,:), allocatable :: grad_fourier_so3_r, grad_fourier_so3_i
      real(dp), allocatable :: t_g_r(:,:), t_g_i(:,:), t_f_r(:,:), t_f_i(:,:), t_g_f_rr(:,:), t_g_f_ii(:,:)
      integer :: alpha

      logical :: my_do_descriptor, my_do_grad_descriptor, do_two_l_plus_one
      integer :: d, i, j, n, a, b, k, l, m, i_pow, i_coeff, l_n_neighbours, n_i, n_descriptors, n_cross, i_species, j_species, ia, jb, i_desc_i, &
         xml_version, sum_l_n_neighbours, i_pair, i_pair_i
      integer, dimension(3) :: shift_ij
      integer, dimension(:), allocatable :: i_desc
      integer, dimension(:,:), allocatable :: rs_index
      real(dp) :: r_ij, arg_bess, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_lp, exp_p, exp_m, f_cut, df_cut, norm_descriptor_i
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:,:), allocatable, save :: radial_fun, radial_coefficient, grad_radial_fun, grad_radial_coefficient, grad_descriptor_i
      real(dp), dimension(:), allocatable, save :: descriptor_i
      real(dp), dimension(:), allocatable :: global_fourier_so3_r_array, global_fourier_so3_i_array
      type(real_2d_array), dimension(:), allocatable :: global_grad_fourier_so3_r_array, global_grad_fourier_so3_i_array
      integer, dimension(116) :: species_map
!$omp threadprivate(radial_fun, radial_coefficient, grad_radial_fun, grad_radial_coefficient)
!$omp threadprivate(fourier_so3_r, fourier_so3_i)
!$omp threadprivate(SphericalY_ij,grad_SphericalY_ij)
!$omp threadprivate(descriptor_i, grad_descriptor_i)

      INIT_ERROR(error)

      call system_timer('soap_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("soap_calc: descriptor object not initialised", error)
      endif

      species_map = 0
      do i_species = 1, this%n_species
         if(this%species_Z(i_species) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i_species)) = i_species
         endif
      enddo

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      has_atom_mask_name = .false.
      atom_mask_pointer => null()
      xml_version = 1423143769

      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
            help_string="Name of a logical property in the atoms object. For atoms where this property is " // &
            "true, descriptors are calculated.")

         call param_register(params, 'xml_version', '1423143769', xml_version, &
            help_string="Version of GAP the XML potential file was created")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='soap_calc args_str')) then
            RAISE_ERROR("soap_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("soap_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      do_two_l_plus_one = (xml_version >= 1423143769)

      allocate(rs_index(2,this%n_max*this%n_species))
      i = 0
      do i_species = 1, this%n_species
         do a = 1, this%n_max
            i = i + 1
            rs_index(:,i) = (/a,i_species/)
         enddo
      enddo

      call finalise(descriptor_out)

      d = soap_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      allocate(i_desc(at%N))

!$omp parallel default(none) shared(this,my_do_grad_descriptor,d) private(i_species, a, l)
      allocate(descriptor_i(d))
      if(my_do_grad_descriptor) allocate(grad_descriptor_i(d,3))

      allocate(radial_fun(0:this%l_max, this%n_max), radial_coefficient(0:this%l_max, this%n_max)) 
      !SPEED allocate(fourier_so3(0:this%l_max,this%n_max,this%n_species), SphericalY_ij(0:this%l_max))
      allocate(fourier_so3_r(0:this%l_max,this%n_max,this%n_species), fourier_so3_i(0:this%l_max,this%n_max,this%n_species), SphericalY_ij(0:this%l_max))

      if(my_do_grad_descriptor) then
         allocate(grad_radial_fun(0:this%l_max, this%n_max), grad_radial_coefficient(0:this%l_max, this%n_max))
         allocate(grad_SphericalY_ij(0:this%l_max))
      endif

      do i_species = 1, this%n_species
         do a = 1, this%n_max
            do l = 0, this%l_max
               !SPEED allocate(fourier_so3(l,a,i_species)%m(-l:l))
               !SPEED fourier_so3(l,a,i_species)%m(:) = CPLX_ZERO
               allocate(fourier_so3_r(l,a,i_species)%m(-l:l))
               allocate(fourier_so3_i(l,a,i_species)%m(-l:l))
               fourier_so3_r(l,a,i_species)%m(:) = 0.0_dp
               fourier_so3_i(l,a,i_species)%m(:) = 0.0_dp
            enddo
         enddo
      enddo

      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
         if(my_do_grad_descriptor) allocate(grad_SphericalY_ij(l)%mm(3,-l:l))
      enddo
!$omp end parallel

      i_desc = 0
      i_desc_i = 0
      do i = 1, at%N

         if( .not. any( at%Z(i) == this%Z ) .and. .not. any(this%Z == 0) ) cycle

         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc_i = i_desc_i + 1
         i_desc(i) = i_desc_i

         if(.not. this%global) then ! atomic SOAP
            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc_i)%data(d))
               !slow, no need
               !descriptor_out%x(i_desc_i)%data = 0.0_dp
               allocate(descriptor_out%x(i_desc_i)%ci(1))
               descriptor_out%x(i_desc_i)%has_data = .false.
               descriptor_out%x(i_desc_i)%covariance_cutoff = 1.0_dp
            endif
            if(my_do_grad_descriptor) then
               l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

               allocate(descriptor_out%x(i_desc_i)%grad_data(d,3,0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%ii(0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%pos(3,0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%has_grad_data(0:l_n_neighbours))
               ! slow, no need
               ! descriptor_out%x(i_desc_i)%grad_data = 0.0_dp
               descriptor_out%x(i_desc_i)%grad_data(:,:,0) = 0.0_dp
               descriptor_out%x(i_desc_i)%ii = 0
               descriptor_out%x(i_desc_i)%pos = 0.0_dp
               descriptor_out%x(i_desc_i)%has_grad_data = .false.

               allocate(descriptor_out%x(i_desc_i)%grad_covariance_cutoff(3,0:l_n_neighbours))
               descriptor_out%x(i_desc_i)%grad_covariance_cutoff = 0.0_dp
            endif
         endif
      enddo

      allocate( &
         global_fourier_so3_r_array((this%l_max+1)**2 * this%n_max * this%n_species), &
         global_fourier_so3_i_array((this%l_max+1)**2 * this%n_max * this%n_species), &
         global_grad_fourier_so3_r_array( count(i_desc/=0) ), &
         global_grad_fourier_so3_i_array( count(i_desc/=0) ) )



      if(this%global) then
         if(my_do_descriptor) then
            allocate(descriptor_out%x(1)%data(d))
            if( any(this%Z == 0) ) then
               allocate(descriptor_out%x(1)%ci(at%N))
               descriptor_out%x(1)%ci(:) = (/ (i, i=1, at%N) /)
            else
               allocate( descriptor_out%x(1)%ci( count( (/(any(at%Z(i)==this%Z),i=1,at%N)/) ) ) )
               forall(i=1:at%N, any(at%Z(i) == this%Z)) descriptor_out%x(1)%ci(i_desc(i)) = i
            endif
            descriptor_out%x(1)%has_data = .true.
            descriptor_out%x(1)%covariance_cutoff = 1.0_dp
         endif ! my_do_descriptor
         if(my_do_grad_descriptor) then
            sum_l_n_neighbours = 0
            do i = 1, at%N

               if(i_desc(i) == 0) then
                  cycle
               else
                  i_desc_i = i_desc(i)
               endif

               l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)
               sum_l_n_neighbours = sum_l_n_neighbours + l_n_neighbours + 1 ! include central atom as well!

               allocate( &
               global_grad_fourier_so3_r_array(i_desc_i)%x(0:this%l_max,this%n_max,l_n_neighbours), &
               global_grad_fourier_so3_i_array(i_desc_i)%x(0:this%l_max,this%n_max,l_n_neighbours) )

               do n_i = 1, l_n_neighbours
                  do a = 1, this%n_max
                     do l = 0, this%l_max
                        allocate( &
                           global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm(3,-l:l), &
                           global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm(3,-l:l) )
                        global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm = 0.0_dp
                        global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm = 0.0_dp
                     enddo ! l
                  enddo ! a
               enddo ! n_i
            enddo ! i

            allocate(descriptor_out%x(1)%grad_data(d,3,sum_l_n_neighbours))
            allocate(descriptor_out%x(1)%ii(sum_l_n_neighbours))
            allocate(descriptor_out%x(1)%pos(3,sum_l_n_neighbours))
            allocate(descriptor_out%x(1)%has_grad_data(sum_l_n_neighbours))

            allocate(descriptor_out%x(1)%grad_covariance_cutoff(3,sum_l_n_neighbours))
            descriptor_out%x(1)%grad_covariance_cutoff = 0.0_dp
         endif ! my_do_grad_descriptor

         global_fourier_so3_r_array = 0.0_dp
         global_fourier_so3_i_array = 0.0_dp
      endif ! this%global 

!$omp parallel do schedule(dynamic) default(none) shared(this, at, descriptor_out, my_do_descriptor, my_do_grad_descriptor, d, i_desc, species_map, rs_index, do_two_l_plus_one) &
!$omp shared(global_grad_fourier_so3_r_array, global_grad_fourier_so3_i_array) &
!$omp private(i, j, i_species, j_species, a, b, l, m, n, n_i, r_ij, u_ij, d_ij, shift_ij, i_pow, i_coeff, ia, jb, alpha, i_desc_i) &
!$omp private(grad_fourier_so3_r,grad_fourier_so3_i,t_g_r, t_g_i, t_f_r, t_f_i, t_g_f_rr, t_g_f_ii) &
!$omp private(f_cut, df_cut, arg_bess, exp_p, exp_m, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lp, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, norm_descriptor_i) &
!$omp reduction(+:global_fourier_so3_r_array,global_fourier_so3_i_array)
      do i = 1, at%N

         if(i_desc(i) == 0) then
            cycle
         else
            i_desc_i = i_desc(i)
         endif

         if(.not.this%global) then
            if(my_do_descriptor) then
               descriptor_out%x(i_desc_i)%ci(1) = i
               descriptor_out%x(i_desc_i)%has_data = .true.
            endif
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc_i)%ii(0) = i
               descriptor_out%x(i_desc_i)%pos(:,0) = at%pos(:,i)
               descriptor_out%x(i_desc_i)%has_grad_data(0) = .true.
            endif
         endif
         if(my_do_grad_descriptor) then
            !SPEED allocate( grad_fourier_so3(0:this%l_max,this%n_max,n_neighbours(at,i,max_dist=this%cutoff)) )
            allocate( grad_fourier_so3_r(0:this%l_max,this%n_max,n_neighbours(at,i,max_dist=this%cutoff)) )
            allocate( grad_fourier_so3_i(0:this%l_max,this%n_max,n_neighbours(at,i,max_dist=this%cutoff)) )
         endif

         !do a = 1, this%n_max
         !   radial_fun(0,a) = exp( -this%alpha * this%r_basis(a)**2 ) !* this%r_basis(a)
         !enddo
         !radial_coefficient(0,:) = matmul( radial_fun(0,:), this%transform_basis )
         radial_fun(0,:) = 0.0_dp
         radial_fun(0,1) = 1.0_dp
         radial_coefficient(0,:) = matmul( radial_fun(0,:), this%cholesky_overlap_basis)

         do i_species = 1, this%n_species
            do a = 1, this%n_max
               !SPEED fourier_so3(0,a,i_species)%m(0) = radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/))
               if( this%central_reference_all_species .or. this%species_Z(i_species) == at%Z(i) .or. this%species_Z(i_species) == 0 ) then
                  fourier_so3_r(0,a,i_species)%m(0) = this%central_weight * real(radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/)), dp)
                  fourier_so3_i(0,a,i_species)%m(0) = this%central_weight * aimag(radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/)))
               else
                  fourier_so3_i(0,a,i_species)%m(0) = 0.0_dp
                  fourier_so3_r(0,a,i_species)%m(0) = 0.0_dp
               endif

               do l = 1, this%l_max
                  !SPEED fourier_so3(l,a,i_species)%m(:) = CPLX_ZERO
                  fourier_so3_r(l,a,i_species)%m(:) = 0.0_dp
                  fourier_so3_i(l,a,i_species)%m(:) = 0.0_dp
               enddo
            enddo
         enddo

! soap_calc 20 takes 0.0052 s
! call system_timer("soap_calc 20")
         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1

            i_species = species_map(at%Z(j))
            if( i_species == 0 ) cycle

            if(.not. this%global .and. my_do_grad_descriptor) then
               descriptor_out%x(i_desc_i)%ii(n_i) = j
               descriptor_out%x(i_desc_i)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc_i)%has_grad_data(n_i) = .true.
            endif

            f_cut = coordination_function(r_ij,this%cutoff, this%cutoff_transition_width)
            if(my_do_grad_descriptor) then
               df_cut = dcoordination_function(r_ij,this%cutoff, this%cutoff_transition_width)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     !SPEED allocate(grad_fourier_so3(l,a,n_i)%mm(3,-l:l))
                     !SPEED grad_fourier_so3(l,a,n_i)%mm(:,:) = CPLX_ZERO
                     allocate(grad_fourier_so3_r(l,a,n_i)%mm(3,-l:l))
                     allocate(grad_fourier_so3_i(l,a,n_i)%mm(3,-l:l))
                     grad_fourier_so3_r(l,a,n_i)%mm(:,:) = 0.0_dp
                     grad_fourier_so3_i(l,a,n_i)%mm(:,:) = 0.0_dp
                  enddo
               enddo
            endif

            do a = 1, this%n_max
               arg_bess = 2.0_dp * this%alpha * r_ij * this%r_basis(a)
               exp_p = exp( -this%alpha*( r_ij + this%r_basis(a) )**2 )
               exp_m = exp( -this%alpha*( r_ij - this%r_basis(a) )**2 )

               do l = 0, this%l_max
                  if( l == 0 ) then
                     if(arg_bess == 0.0_dp) then
                        !mo_spher_bess_fi_ki_l = 1.0_dp
                        mo_spher_bess_fi_ki_l = exp( -this%alpha * (this%r_basis(a)**2 + r_ij**2) )
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp
                     else
                        !mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                        !mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                        mo_spher_bess_fi_ki_lm = 0.5_dp * (exp_m + exp_p) / arg_bess
                        mo_spher_bess_fi_ki_l  = 0.5_dp * (exp_m - exp_p) / arg_bess
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                     endif
                  else
                     if(arg_bess == 0.0_dp) then
                        mo_spher_bess_fi_ki_l = 0.0_dp
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp
                     else
                        mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                        mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                        if(my_do_grad_descriptor) then
                           mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                           mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                        else
                           mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess
                        endif
                     endif
                  endif

                  !radial_fun(l,a) = exp( -this%alpha * (this%r_basis(a)**2 + r_ij**2) ) * mo_spher_bess_fi_ki_l !* this%r_basis(a)
                  radial_fun(l,a) = mo_spher_bess_fi_ki_l !* this%r_basis(a)
                  if(my_do_grad_descriptor) grad_radial_fun(l,a) = -2.0_dp * this%alpha * r_ij * mo_spher_bess_fi_ki_l + &
                     l*mo_spher_bess_fi_ki_l / r_ij + mo_spher_bess_fi_ki_lp * 2.0_dp * this%alpha * this%r_basis(a)

               enddo
            enddo

            radial_coefficient = matmul( radial_fun, this%transform_basis )
            if(my_do_grad_descriptor) grad_radial_coefficient = matmul( grad_radial_fun, this%transform_basis ) * f_cut + radial_coefficient * df_cut
            radial_coefficient = radial_coefficient * f_cut

            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian(l,m,d_ij)
                  if(my_do_grad_descriptor) grad_SphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian(l,m,d_ij)
               enddo
            enddo

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = -l, l
                     !SPEED fourier_so3(l,a,i_species)%m(m) = fourier_so3(l,a,i_species)%m(m) + radial_coefficient(l,a) * SphericalY_ij(l)%m(m)
                     !SPEED if(my_do_grad_descriptor) grad_fourier_so3(l,a,n_i)%mm(:,m) = grad_fourier_so3(l,a,n_i)%mm(:,m) + &
                     !SPEED    grad_radial_coefficient(l,a) * SphericalY_ij(l)%m(m) * u_ij + radial_coefficient(l,a) * grad_SphericalY_ij(l)%mm(:,m)
                     fourier_so3_r(l,a,i_species)%m(m) = fourier_so3_r(l,a,i_species)%m(m) + real(radial_coefficient(l,a) * SphericalY_ij(l)%m(m), dp)
                     fourier_so3_i(l,a,i_species)%m(m) = fourier_so3_i(l,a,i_species)%m(m) + aimag(radial_coefficient(l,a) * SphericalY_ij(l)%m(m))
                     if(my_do_grad_descriptor) then
                        grad_fourier_so3_r(l,a,n_i)%mm(:,m) = grad_fourier_so3_r(l,a,n_i)%mm(:,m) + &
                           real(grad_radial_coefficient(l,a) * SphericalY_ij(l)%m(m) * u_ij + radial_coefficient(l,a) * grad_SphericalY_ij(l)%mm(:,m), dp)
                        grad_fourier_so3_i(l,a,n_i)%mm(:,m) = grad_fourier_so3_i(l,a,n_i)%mm(:,m) + &
                           aimag(grad_radial_coefficient(l,a) * SphericalY_ij(l)%m(m) * u_ij + radial_coefficient(l,a) * grad_SphericalY_ij(l)%mm(:,m))
                     endif ! my_do_grad_descriptor
                  enddo ! m
               enddo ! l
            enddo ! a

         enddo ! n
! call system_timer("soap_calc 20")

         if(this%global .and. my_do_grad_descriptor) then
            global_grad_fourier_so3_r_array(i_desc_i)%x = grad_fourier_so3_r
            global_grad_fourier_so3_i_array(i_desc_i)%x = grad_fourier_so3_i
            !do n_i = lbound(grad_fourier_so3_r,3), ubound(grad_fourier_so3_r,3)
            !   do a = lbound(grad_fourier_so3_r,2), ubound(grad_fourier_so3_r,2)
            !      do l = lbound(grad_fourier_so3_r,1), ubound(grad_fourier_so3_r,1)
            !         global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm = grad_fourier_so3_r(l,a,n_i)%mm
            !         global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm = grad_fourier_so3_i(l,a,n_i)%mm
            !      enddo ! l
            !   enddo ! a
            !enddo ! n_i
         endif

         if(this%global) then
            i_coeff = 0
            do ia = 1, this%n_species*this%n_max
               a = rs_index(1,ia)
               i_species = rs_index(2,ia)
               do l = 0, this%l_max
                  global_fourier_so3_r_array(i_coeff+1:i_coeff+2*l+1) = global_fourier_so3_r_array(i_coeff+1:i_coeff+2*l+1) + fourier_so3_r(l,a,i_species)%m(:)
                  global_fourier_so3_i_array(i_coeff+1:i_coeff+2*l+1) = global_fourier_so3_i_array(i_coeff+1:i_coeff+2*l+1) + fourier_so3_i(l,a,i_species)%m(:)
                  i_coeff = i_coeff + 2*l+1
               enddo
            enddo
         endif

         i_pow = 0
         do ia = 1, this%n_species*this%n_max
            a = rs_index(1,ia)
            i_species = rs_index(2,ia)
            do jb = 1, ia
               b = rs_index(1,jb)
               j_species = rs_index(2,jb)

               if(this%diagonal_radial .and. a /= b) cycle

               do l = 0, this%l_max
                  i_pow = i_pow + 1
                  !SPEED descriptor_i(i_pow) = real( dot_product(fourier_so3(l,a,i_species)%m, fourier_so3(l,b,j_species)%m) )
                  descriptor_i(i_pow) = dot_product(fourier_so3_r(l,a,i_species)%m, fourier_so3_r(l,b,j_species)%m) + dot_product(fourier_so3_i(l,a,i_species)%m, fourier_so3_i(l,b,j_species)%m)
                  if(do_two_l_plus_one) descriptor_i(i_pow) = descriptor_i(i_pow) / sqrt(2.0_dp * l + 1.0_dp)
                  if( ia /= jb ) descriptor_i(i_pow) = descriptor_i(i_pow) * SQRT_TWO
               enddo !l
            enddo !jb
         enddo !ia

         descriptor_i(d) = 0.0_dp
         norm_descriptor_i = sqrt(dot_product(descriptor_i,descriptor_i))

         if(.not. this%global .and. my_do_descriptor) then
            if(this%normalise) then
               descriptor_out%x(i_desc_i)%data = descriptor_i / norm_descriptor_i
            else
               descriptor_out%x(i_desc_i)%data = descriptor_i
            endif

            descriptor_out%x(i_desc_i)%data(d) = this%covariance_sigma0
         endif

         if(my_do_grad_descriptor) then
! soap_calc 33 takes 0.047 s
! call system_timer("soap_calc 33")
	    allocate(t_g_r(this%n_max*3, 2*this%l_max+1), t_g_i(this%n_max*3, 2*this%l_max+1))
	    allocate(t_f_r(this%n_max*this%n_species, 2*this%l_max+1), t_f_i(this%n_max*this%n_species, 2*this%l_max+1))
	    allocate(t_g_f_rr(this%n_max*3, this%n_max*this%n_species), t_g_f_ii(this%n_max*3, this%n_max*this%n_species))
            !do n_i = 1, n_neighbours(at,i,max_dist=this%cutoff)

            n_i = 0
            do n = 1, n_neighbours(at,i)
               j = neighbour(at, i, n, distance = r_ij)
               if( r_ij >= this%cutoff ) cycle

               n_i = n_i + 1

               if( species_map(at%Z(j)) == 0 ) cycle

               i_pow = 0
               grad_descriptor_i = 0.0_dp

               !SPEED do ia = 1, this%n_species*this%n_max
               !SPEED    a = rs_index(1,ia) 
               !SPEED    i_species = rs_index(2,ia)
               !SPEED    do jb = 1, ia
               !SPEED       b = rs_index(1,jb)
               !SPEED       j_species = rs_index(2,jb)
               !SPEED       do l = 0, this%l_max
               !SPEED          i_pow = i_pow + 1
               !SPEED          if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) grad_descriptor_i(i_pow,:) = grad_descriptor_i(i_pow,:) + real( matmul(conjg(grad_fourier_so3(l,a,n_i)%mm),fourier_so3(l,b,j_species)%m) )
               !SPEED          if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) grad_descriptor_i(i_pow,:) = grad_descriptor_i(i_pow,:) + real( matmul(grad_fourier_so3(l,b,n_i)%mm,conjg(fourier_so3(l,a,i_species)%m)) )
               !SPEED          !grad_descriptor_i(i_pow,:) = real( matmul(conjg(grad_fourier_so3(l,a,n_i)%mm),fourier_so3(l,b)%m) + matmul(grad_fourier_so3(l,b,n_i)%mm,conjg(fourier_so3(l,a)%m)) )
               !SPEED          if( ia /= jb ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
               !SPEED       enddo !l
               !SPEED    enddo !jb
               !SPEED enddo !ia

               !SPEED do ia = 1, this%n_species*this%n_max
               !SPEED    a = rs_index(1,ia) 
               !SPEED    i_species = rs_index(2,ia)
               !SPEED    do jb = 1, ia
               !SPEED       b = rs_index(1,jb)
               !SPEED       j_species = rs_index(2,jb)
               !SPEED       do l = 0, this%l_max
               !SPEED          i_pow = i_pow + 1
               !SPEED          if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) grad_descriptor_i(i_pow,:) = grad_descriptor_i(i_pow,:) + &
               !SPEED             matmul(grad_fourier_so3_r(l,a,n_i)%mm,fourier_so3_r(l,b,j_species)%m) + matmul(grad_fourier_so3_i(l,a,n_i)%mm,fourier_so3_i(l,b,j_species)%m)
               !SPEED          if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) grad_descriptor_i(i_pow,:) = grad_descriptor_i(i_pow,:) + &
               !SPEED             matmul(grad_fourier_so3_r(l,b,n_i)%mm,fourier_so3_r(l,a,i_species)%m) + matmul(grad_fourier_so3_i(l,b,n_i)%mm,fourier_so3_i(l,a,i_species)%m)
               !SPEED          if( ia /= jb ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
               !SPEED       enddo !l
               !SPEED    enddo !jb
               !SPEED enddo !ia

               do l=0, this%l_max
                  do a = 1, this%n_max
                     do alpha=1, 3
                	t_g_r(3*(a-1)+alpha, 1:2*l+1) = grad_fourier_so3_r(l,a,n_i)%mm(alpha,-l:l)
                	t_g_i(3*(a-1)+alpha, 1:2*l+1) = grad_fourier_so3_i(l,a,n_i)%mm(alpha,-l:l)
                     enddo
                  enddo
                  do ia = 1, this%n_species*this%n_max
                     a = rs_index(1,ia)
                     i_species = rs_index(2,ia)
                     
                     t_f_r(ia, 1:2*l+1) = fourier_so3_r(l,a,i_species)%m(-l:l)
                     t_f_i(ia, 1:2*l+1) = fourier_so3_i(l,a,i_species)%m(-l:l)
                  enddo
                  call dgemm('N','T',this%n_max*3, this%n_max*this%n_species, 2*l+1, 1.0_dp, &
                     t_g_r(1,1), size(t_g_r,1), t_f_r(1,1), size(t_f_r,1), 0.0_dp, t_g_f_rr(1,1), size(t_g_f_rr, 1))
                  call dgemm('N','T',this%n_max*3, this%n_max*this%n_species, 2*l+1, 1.0_dp, &
                     t_g_i(1,1), size(t_g_i,1), t_f_i(1,1), size(t_f_i,1), 0.0_dp, t_g_f_ii(1,1), size(t_g_f_ii, 1))
                  !t_g_f_rr = matmul(t_g_r,transpose(t_f_r))
                  !t_g_f_ii = matmul(t_g_i,transpose(t_f_i))
               
                  i_pow = l+1
                  do ia = 1, this%n_species*this%n_max
                     a = rs_index(1,ia)
                     i_species = rs_index(2,ia)
                     do jb = 1, ia !this%n_species*this%n_max !ia
                        b = rs_index(1,jb)
                        j_species = rs_index(2,jb)

                        if(this%diagonal_radial .and. a /= b) cycle
               
                        if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) + t_g_f_rr(3*(a-1)+1:3*a,jb) + t_g_f_ii(3*(a-1)+1:3*a,jb)
                        if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) + t_g_f_rr(3*(b-1)+1:3*b,ia) + t_g_f_ii(3*(b-1)+1:3*b,ia)
               

                        if(do_two_l_plus_one) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) / sqrt(2.0_dp * l + 1.0_dp)
                        if( ia /= jb ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
                        i_pow = i_pow + this%l_max+1
                     enddo
                  enddo

                  !do a = 1, this%n_max
                  !   do b = 1, a
                  !      grad_descriptor_i(i_pow, 1:3) = t_g_f_rr(3*(a-1)+1:3*a,b) + t_g_f_ii(3*(a-1)+1:3*a,b) + &
                  !                                      t_g_f_rr(3*(b-1)+1:3*b,a) + t_g_f_ii(3*(b-1)+1:3*b,a)
                  !      if( a /= b ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
                  !      i_pow = i_pow + this%l_max+1
                  !   end do
                  !end do

               end do !l

               grad_descriptor_i(d, 1:3) = 0.0_dp

               if(.not. this%global) then
                  if( this%normalise ) then
                     descriptor_out%x(i_desc_i)%grad_data(:,:,n_i) = grad_descriptor_i / norm_descriptor_i
                     do k = 1, 3
                        descriptor_out%x(i_desc_i)%grad_data(:,k,n_i) = descriptor_out%x(i_desc_i)%grad_data(:,k,n_i) - descriptor_i * dot_product(descriptor_i,grad_descriptor_i(:,k)) / norm_descriptor_i**3
                     enddo
                  else
                     descriptor_out%x(i_desc_i)%grad_data(:,:,n_i) = grad_descriptor_i
                  endif

                  descriptor_out%x(i_desc_i)%grad_data(:,:,0) = descriptor_out%x(i_desc_i)%grad_data(:,:,0) - descriptor_out%x(i_desc_i)%grad_data(:,:,n_i)
               endif
            enddo !ni
	    deallocate(t_f_r, t_f_i)
	    deallocate(t_g_r, t_g_i)
	    deallocate(t_g_f_rr, t_g_f_ii)
! call system_timer("soap_calc 33")

            do n_i = 1, n_neighbours(at,i,max_dist=this%cutoff)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     !SPEED deallocate(grad_fourier_so3(l,a,n_i)%mm)
                     if(allocated(grad_fourier_so3_r(l,a,n_i)%mm)) deallocate(grad_fourier_so3_r(l,a,n_i)%mm)
                     if(allocated(grad_fourier_so3_i(l,a,n_i)%mm)) deallocate(grad_fourier_so3_i(l,a,n_i)%mm)
                  enddo
               enddo
            enddo
            !SPEED deallocate(grad_fourier_so3)
            deallocate(grad_fourier_so3_r)
            deallocate(grad_fourier_so3_i)
         endif

      enddo ! i
!$omp end parallel do

      !SPEED if(allocated(fourier_so3)) then
      !SPEED    do i_species = 1, this%n_species
      !SPEED       do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
      !SPEED          do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
      !SPEED             deallocate(fourier_so3(l,a,i_species)%m)
      !SPEED          enddo
      !SPEED       enddo
      !SPEED    enddo
      !SPEED    deallocate(fourier_so3)
      !SPEED endif

!$omp parallel default(none) private(i_species, a, l)
      if(allocated(fourier_so3_r)) then
         do i_species = lbound(fourier_so3_r,3), ubound(fourier_so3_r,3)
            do a = lbound(fourier_so3_r,2), ubound(fourier_so3_r,2)
               do l = lbound(fourier_so3_r,1), ubound(fourier_so3_r,1)
                  deallocate(fourier_so3_r(l,a,i_species)%m)
               enddo
            enddo
         enddo
         deallocate(fourier_so3_r)
      endif
      if(allocated(fourier_so3_i)) then
         do i_species = lbound(fourier_so3_i,3), ubound(fourier_so3_i,3)
            do a = lbound(fourier_so3_i,2), ubound(fourier_so3_i,2)
               do l = lbound(fourier_so3_i,1), ubound(fourier_so3_i,1)
                  deallocate(fourier_so3_i(l,a,i_species)%m)
               enddo
            enddo
         enddo
         deallocate(fourier_so3_i)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(grad_SphericalY_ij)) then
         do l = lbound(grad_SphericalY_ij,1), ubound(grad_SphericalY_ij,1)
            deallocate(grad_SphericalY_ij(l)%mm)
         enddo
         deallocate(grad_SphericalY_ij)
      endif

      if(allocated(radial_fun)) deallocate(radial_fun)
      if(allocated(radial_coefficient)) deallocate(radial_coefficient)
      if(allocated(grad_radial_fun)) deallocate(grad_radial_fun)
      if(allocated(grad_radial_coefficient)) deallocate(grad_radial_coefficient)
      if(allocated(descriptor_i)) deallocate(descriptor_i)
      if(allocated(grad_descriptor_i)) deallocate(grad_descriptor_i)
!$omp end parallel

      if(this%global) then
         allocate(global_fourier_so3_r(0:this%l_max,this%n_max,this%n_species), global_fourier_so3_i(0:this%l_max,this%n_max,this%n_species), &
            descriptor_i(d) )

         i_coeff = 0
         do ia = 1, this%n_species*this%n_max
            a = rs_index(1,ia)
            i_species = rs_index(2,ia)
            do l = 0, this%l_max
               allocate(global_fourier_so3_r(l,a,i_species)%m(-l:l))
               allocate(global_fourier_so3_i(l,a,i_species)%m(-l:l))
               global_fourier_so3_r(l,a,i_species)%m(:) = global_fourier_so3_r_array(i_coeff+1:i_coeff+2*l+1)
               global_fourier_so3_i(l,a,i_species)%m(:) = global_fourier_so3_i_array(i_coeff+1:i_coeff+2*l+1)
               i_coeff = i_coeff + 2*l+1
            enddo
         enddo

         i_pow = 0
         do ia = 1, this%n_species*this%n_max
            a = rs_index(1,ia)
            i_species = rs_index(2,ia)
            do jb = 1, ia
               b = rs_index(1,jb)
               j_species = rs_index(2,jb)

               if(this%diagonal_radial .and. a /= b) cycle

               do l = 0, this%l_max
                  i_pow = i_pow + 1
                  descriptor_i(i_pow) = &
                     dot_product(global_fourier_so3_r(l,a,i_species)%m, global_fourier_so3_r(l,b,j_species)%m) + &
                     dot_product(global_fourier_so3_i(l,a,i_species)%m, global_fourier_so3_i(l,b,j_species)%m)

                  !if( ia /= jb ) descriptor_out%global_data(i_pow) = descriptor_out%global_data(i_pow) * SQRT_TWO
                  if(do_two_l_plus_one) descriptor_i(i_pow) = descriptor_i(i_pow) / sqrt(2.0_dp * l + 1.0_dp)
                  if( ia /= jb ) descriptor_i(i_pow) = descriptor_i(i_pow) * SQRT_TWO
               enddo !l

            enddo !jb
         enddo !ia
         descriptor_i(d) = 0.0_dp
         norm_descriptor_i = sqrt(dot_product(descriptor_i,descriptor_i))
         if(my_do_descriptor) then
            if(this%normalise) then
               descriptor_out%x(1)%data = descriptor_i / norm_descriptor_i
            else
               descriptor_out%x(1)%data = descriptor_i
            endif
            descriptor_out%x(1)%data(d) = this%covariance_sigma0
         endif

         if(my_do_grad_descriptor) then
	    allocate(t_g_r(this%n_max*3, 2*this%l_max+1), t_g_i(this%n_max*3, 2*this%l_max+1))
	    allocate(t_f_r(this%n_max*this%n_species, 2*this%l_max+1), t_f_i(this%n_max*this%n_species, 2*this%l_max+1))
	    allocate(t_g_f_rr(this%n_max*3, this%n_max*this%n_species), t_g_f_ii(this%n_max*3, this%n_max*this%n_species))
            allocate(grad_descriptor_i(d,3))

            i_pair = 0
            do i = 1, at%N

               if(i_desc(i) == 0) then
                  cycle
               else
                  i_desc_i = i_desc(i)
               endif

               i_pair = i_pair + 1
               i_pair_i = i_pair ! accumulates \frac{ \partial p^{(j)} }{ \partial r_{ji\alpha} }

               descriptor_out%x(1)%ii(i_pair_i) = i
               descriptor_out%x(1)%pos(:,i_pair_i) = 0.0_dp
               descriptor_out%x(1)%has_grad_data(i_pair_i) = .true.
               descriptor_out%x(1)%grad_data(:,:,i_pair_i) = 0.0_dp

               n_i = 0
               do n = 1, n_neighbours(at,i)
                  j = neighbour(at, i, n, distance = r_ij, diff = d_ij)
                  if( r_ij >= this%cutoff ) cycle

                  n_i = n_i + 1
                  i_pair = i_pair + 1 ! \frac{ \partial p^{(i)} }{ \partial r_{ij\alpha} }

                  descriptor_out%x(1)%ii(i_pair) = j
                  descriptor_out%x(1)%pos(:,i_pair) = d_ij
                  descriptor_out%x(1)%has_grad_data(i_pair) = .true.

                  i_pow = 0
                  grad_descriptor_i = 0.0_dp

                  do l=0, this%l_max
                     do a = 1, this%n_max
                        do alpha=1, 3
                           t_g_r(3*(a-1)+alpha, 1:2*l+1) = global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm(alpha,-l:l)
                           t_g_i(3*(a-1)+alpha, 1:2*l+1) = global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm(alpha,-l:l)
                        enddo
                     enddo
                     do ia = 1, this%n_species*this%n_max
                        a = rs_index(1,ia)
                        i_species = rs_index(2,ia)
                        
                        t_f_r(ia, 1:2*l+1) = global_fourier_so3_r(l,a,i_species)%m(-l:l)
                        t_f_i(ia, 1:2*l+1) = global_fourier_so3_i(l,a,i_species)%m(-l:l)
                     enddo
                     call dgemm('N','T',this%n_max*3, this%n_max*this%n_species, 2*l+1, 1.0_dp, &
                        t_g_r(1,1), size(t_g_r,1), t_f_r(1,1), size(t_f_r,1), 0.0_dp, t_g_f_rr(1,1), size(t_g_f_rr, 1))
                     call dgemm('N','T',this%n_max*3, this%n_max*this%n_species, 2*l+1, 1.0_dp, &
                        t_g_i(1,1), size(t_g_i,1), t_f_i(1,1), size(t_f_i,1), 0.0_dp, t_g_f_ii(1,1), size(t_g_f_ii, 1))
                     !t_g_f_rr = matmul(t_g_r,transpose(t_f_r))
                     !t_g_f_ii = matmul(t_g_i,transpose(t_f_i))
                  
                     i_pow = l+1
                     do ia = 1, this%n_species*this%n_max
                        a = rs_index(1,ia)
                        i_species = rs_index(2,ia)
                        do jb = 1, ia !this%n_species*this%n_max !ia
                           b = rs_index(1,jb)
                           j_species = rs_index(2,jb)
                  
                           if(this%diagonal_radial .and. a /= b) cycle

                           if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) + t_g_f_rr(3*(a-1)+1:3*a,jb) + t_g_f_ii(3*(a-1)+1:3*a,jb)
                           if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) + t_g_f_rr(3*(b-1)+1:3*b,ia) + t_g_f_ii(3*(b-1)+1:3*b,ia)
                  

                           if(do_two_l_plus_one) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) / sqrt(2.0_dp * l + 1.0_dp)
                           if( ia /= jb ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
                           i_pow = i_pow + this%l_max+1
                        enddo
                     enddo

                  end do !l

                  grad_descriptor_i(d, 1:3) = 0.0_dp
                  if( this%normalise ) then
                     descriptor_out%x(1)%grad_data(:,:,i_pair) = grad_descriptor_i / norm_descriptor_i
                     do k = 1, 3
                        descriptor_out%x(1)%grad_data(:,k,i_pair) = descriptor_out%x(1)%grad_data(:,k,i_pair) - descriptor_i * dot_product(descriptor_i,grad_descriptor_i(:,k)) / norm_descriptor_i**3
                     enddo
                  else
                     descriptor_out%x(1)%grad_data(:,:,i_pair) = grad_descriptor_i
                  endif

                  descriptor_out%x(1)%grad_data(:,:,i_pair_i) = descriptor_out%x(1)%grad_data(:,:,i_pair_i) - descriptor_out%x(1)%grad_data(:,:,i_pair)
               enddo ! n/n_i

            enddo ! i

            deallocate(grad_descriptor_i)
            deallocate(t_f_r, t_f_i)
            deallocate(t_g_r, t_g_i)
            deallocate(t_g_f_rr, t_g_f_ii)
         endif ! my_do_grad_descriptor

         do i_species = lbound(global_fourier_so3_r,3), ubound(global_fourier_so3_r,3)
            do a = lbound(global_fourier_so3_r,2), ubound(global_fourier_so3_r,2)
               do l = lbound(global_fourier_so3_r,1), ubound(global_fourier_so3_r,1)
                  deallocate(global_fourier_so3_r(l,a,i_species)%m)
               enddo
            enddo
         enddo
         deallocate(global_fourier_so3_r)
         do i_species = lbound(global_fourier_so3_i,3), ubound(global_fourier_so3_i,3)
            do a = lbound(global_fourier_so3_i,2), ubound(global_fourier_so3_i,2)
               do l = lbound(global_fourier_so3_i,1), ubound(global_fourier_so3_i,1)
                  deallocate(global_fourier_so3_i(l,a,i_species)%m)
               enddo
            enddo
         enddo
         deallocate(global_fourier_so3_i)

         if(allocated(descriptor_i)) deallocate(descriptor_i)
      endif ! this%global

      if(allocated(global_fourier_so3_r_array)) deallocate(global_fourier_so3_r_array)
      if(allocated(global_fourier_so3_i_array)) deallocate(global_fourier_so3_i_array)

      if(allocated(global_grad_fourier_so3_r_array)) then
         do i_desc_i = lbound(global_grad_fourier_so3_r_array,1), ubound(global_grad_fourier_so3_r_array,1)
            if(allocated(global_grad_fourier_so3_r_array(i_desc_i)%x)) then
               do n_i = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,3), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,3)
                  do a = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,2), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,2)
                     do l = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,1), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,1)
                        if(allocated(global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm)) deallocate(global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm)
                     enddo ! l
                  enddo ! a
               enddo ! n_i
               deallocate(global_grad_fourier_so3_r_array(i_desc_i)%x)
            endif
         enddo ! i_desc_i
         deallocate(global_grad_fourier_so3_r_array)
      endif
      if(allocated(global_grad_fourier_so3_i_array)) then
         do i_desc_i = lbound(global_grad_fourier_so3_i_array,1), ubound(global_grad_fourier_so3_i_array,1)
            if(allocated(global_grad_fourier_so3_i_array(i_desc_i)%x)) then
               do n_i = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,3), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,3)
                  do a = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,2), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,2)
                     do l = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,1), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,1)
                        if(allocated(global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm)) deallocate(global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm)
                     enddo ! l
                  enddo ! a
               enddo ! n_i
               deallocate(global_grad_fourier_so3_i_array(i_desc_i)%x)
            endif
         enddo ! i_desc_i
         deallocate(global_grad_fourier_so3_i_array)
      endif

      if(allocated(rs_index)) deallocate(rs_index)
      if(allocated(i_desc)) deallocate(i_desc)

      call system_timer('soap_calc')

   endsubroutine soap_calc

   subroutine AN_monomer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(AN_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m
      integer, dimension(3) :: shift_ij, shift_ik
      real(dp) :: r_ij, r_ik, r_jk
      real(dp), dimension(3) :: d_ij, d_ik, d_jk, u_ij, u_jk

      INIT_ERROR(error)

      call system_timer('AN_monomer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("AN_monomer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='AN_monomer_calc args_str')) then
            RAISE_ERROR("AN_monomer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("AN_monomer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("AN_monomer_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = AN_monomer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            if(this%do_atomic) then
               allocate(descriptor_out%x(i)%ci(1))
            else
               allocate(descriptor_out%x(i)%ci(this%N))
            endif
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:this%N-1))
            allocate(descriptor_out%x(i)%ii(0:this%N-1))
            allocate(descriptor_out%x(i)%pos(3,0:this%N-1))
            allocate(descriptor_out%x(i)%has_grad_data(0:this%N-1))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,0:this%N-1))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      if(at%N /= this%N) then
         RAISE_ERROR("AN_monomer_calc: number of atoms is "//at%N//" instead of "//this%N,error)
      endif

      do i = 1, at%N

         i_desc = 0

         if(my_do_descriptor) then
            if(this%do_atomic) then
               descriptor_out%x(i)%ci(1) = i
            else
               descriptor_out%x(i)%ci(:) = (/(m,m=1,this%N)/)
            endif
         endif

         if(my_do_grad_descriptor) then
            descriptor_out%x(i)%ii(0) = i
            descriptor_out%x(i)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i)%has_grad_data(:) = .true.
         endif

         do n = 1, n_neighbours(at,i)
            j = neighbour(at,i,n,distance=r_ij, cosines=u_ij, shift=shift_ij)

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               descriptor_out%x(i)%has_data = .true.
               descriptor_out%x(i)%data(i_desc) = r_ij
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i)%ii(n) = j
               descriptor_out%x(i)%pos(:,n) = at%pos(:,j) + matmul(at%lattice,shift_ij)

               descriptor_out%x(i)%grad_data(i_desc,:,n) =  u_ij
               descriptor_out%x(i)%grad_data(i_desc,:,0) = -u_ij
            endif

            do m = 1, n_neighbours(at,i)
               if(n >= m) cycle

               k = neighbour(at,i,m,distance=r_ik, shift=shift_ik)

               d_jk = ( at%pos(:,j) + matmul(at%lattice,shift_ij) ) - ( at%pos(:,k) + matmul(at%lattice,shift_ik) )
               r_jk = norm(d_jk)
               u_jk = d_jk / r_jk

               i_desc = i_desc + 1
               if(my_do_descriptor) then
                  descriptor_out%x(i)%has_data = .true.
                  descriptor_out%x(i)%data(i_desc) = r_jk
               endif

               if(my_do_grad_descriptor) then
                  descriptor_out%x(i)%grad_data(i_desc,:,n) =  u_jk 
                  descriptor_out%x(i)%grad_data(i_desc,:,m) = -u_jk 
               endif

            enddo
         enddo

         if(.not. this%do_atomic) exit

      enddo

      call system_timer('AN_monomer_calc')

   endsubroutine AN_monomer_calc

   subroutine general_monomer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(general_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, monomer_size, i, i_atomic, j_atomic, k, start, finish
      integer, dimension(3) :: temp_shift
      real(dp), dimension(:), allocatable :: dist_vec
      real(dp), dimension(:,:), allocatable :: interatomic_distances
      real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
      integer, dimension(:), allocatable :: atomic_index
      integer, dimension(:,:), allocatable :: monomer_index, shifts
      logical, dimension(:), allocatable :: associated_to_monomer


      INIT_ERROR(error)

      call system_timer('general_monomer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("general_monomer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='general_monomer_calc args_str')) then
            RAISE_ERROR("general_monomer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("general_monomer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      monomer_size=size(this%signature)
      d = general_monomer_dimensions(this,error)

      allocate(shifts(monomer_size,3))
      allocate(dist_vec(d))
      allocate(atomic_index(monomer_size))
      allocate(associated_to_monomer(at%N))
      allocate(interatomic_vectors(monomer_size,monomer_size,3))
      allocate(interatomic_distances(monomer_size,monomer_size))
      interatomic_vectors = 0.0_dp
      interatomic_distances = 0.0_dp
      associated_to_monomer=.False.

      call find_general_monomer(at,monomer_index,this%signature,associated_to_monomer,this%cutoff,this%atom_ordercheck,error)
      if(.not. all(associated_to_monomer)) then
         !RAISE_ERROR("general_monomer_calc: not all atoms assigned to a monomer", error)
         call print("Not all atoms can be assigned to a monomer with atomic numbers "//this%signature)
      endif
      n_descriptors = size(monomer_index,2)
      call print("found "//n_descriptors//" monomers", PRINT_VERBOSE)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            allocate(descriptor_out%x(i)%ci(monomer_size))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,monomer_size))
            allocate(descriptor_out%x(i)%ii(monomer_size))
            allocate(descriptor_out%x(i)%pos(3,monomer_size))
            allocate(descriptor_out%x(i)%has_grad_data(monomer_size))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,monomer_size))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo


      do i = 1, n_descriptors

         atomic_index = monomer_index(:,i) !stores the indices of atoms in this monomer
!write(*,*) "THE ATOMS IN THE MONOMER ARE : "// atomic_index

         if(associated(atom_mask_pointer)) then
            if(.not. any(atom_mask_pointer(atomic_index))) then
               cycle
            else
               if(.not. all(atom_mask_pointer(atomic_index))) then
                  RAISE_ERROR("general_monomer_calc: atom mask has to encompass either all or none of the atoms of a monomer",error)
               endif
            endif
         endif

         !calc all positions relative to atom 1
         do i_atomic=2,monomer_size
           temp_shift=0
           interatomic_vectors(1,i_atomic,:) = diff_min_image(at,atomic_index(1),atomic_index(i_atomic),shift=temp_shift)
           shifts(i_atomic,:) = temp_shift
         end do

         !find other relative positions through vector addition
         do j_atomic=2,monomer_size
           do i_atomic=2,j_atomic-1
             interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:)
           end do
         end do

         !Now convert vectors to scalar distances
         do i_atomic=1,monomer_size
           do j_atomic=i_atomic+1,monomer_size
             interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
           end do
         end do
!!$do i_atomic=1,size(interatomic_distances,1)
!!$  write(*,'(6F12.8)') interatomic_distances(i_atomic,:)
!!$end do
         !and convert this NxN matrix into the required vector length N(N-1)/2
         start = 1
         finish = monomer_size-1
         do i_atomic=1,monomer_size-1
           dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:monomer_size)
           start = finish+1
           finish=finish + monomer_size-i_atomic-1
         end do

         if(my_do_descriptor) then
            descriptor_out%x(i)%ci(:) = atomic_index
            descriptor_out%x(i)%has_data = .true.
            descriptor_out%x(i)%data = dist_vec
         endif
         call print("distances: "//dist_vec, PRINT_VERBOSE)
         if(my_do_grad_descriptor) then
!!$write(*,*) "doing grad descriptor"
            descriptor_out%x(i)%ii(:) = atomic_index
!!$do i_atomic=1,at%N
!!$write(*,*) at%pos(:,atomic_index(i_atomic))
!!$end do
            descriptor_out%x(i)%pos(:,1) = at%pos(:,atomic_index(1))
            do i_atomic =2,monomer_size
              descriptor_out%x(i)%pos(:,i_atomic) = at%pos(:,atomic_index(i_atomic)) + matmul(at%lattice,shifts(i_atomic,:))
            end do

            !build the grad_data matrix
            descriptor_out%x(i)%has_grad_data(:) = .true.
            do k=1,d
             !find the pair of atoms contributing to this descriptor
             do i_atomic=1,monomer_size
               do j_atomic=i_atomic+1,monomer_size
                 if (interatomic_distances(i_atomic,j_atomic)==dist_vec(k)) then
                   descriptor_out%x(i)%grad_data(k,:,i_atomic) = -interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! kth descriptor wrt atom i_atomic
                   descriptor_out%x(i)%grad_data(k,:,j_atomic) = -descriptor_out%x(i)%grad_data(k,:,i_atomic)        ! kth descriptor wrt j_atomic
!write(*,*) "descriptor dimension "//k//" wrt atoms "//atomic_index(i_atomic)//" and "//atomic_index(j_atomic)
                 end if
               end do
             end do
            end do

       
         endif

      enddo

      deallocate(shifts)
      deallocate(dist_vec)
      deallocate(atomic_index)
      deallocate(associated_to_monomer)
      deallocate(interatomic_vectors)
      deallocate(interatomic_distances)
      deallocate(monomer_index)
      call system_timer('general_monomer_calc')

   endsubroutine general_monomer_calc

   subroutine com_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)

      type(com_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, use_smooth_cutoff, double_count
      integer :: n_descriptors, dimer_size, i, j, i_atomic, j_atomic, i_desc
      integer :: monomer_one_size, monomer_two_size, n_monomer_one, n_monomer_two, this_pair, diff_loc
      integer, dimension(1) :: unit_array
      real(dp), dimension(3) :: diff_one_two, com_pos_one, com_pos_two, transdirvec
      real(dp) :: dist, primitive_cutoff, primitive_cutoff_grad
      real(dp), dimension(:,:), allocatable :: diffs_one, diffs_two, com_pos_diffs
      real(dp), dimension(:), allocatable :: weight_one, weight_two
      integer, dimension(:), allocatable :: atomic_index, atomic_index_one, atomic_index_two, pairs_diffs_map
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index,  monomer_pairs
      logical, dimension(:), allocatable :: associated_to_monomer


      INIT_ERROR(error)
      use_smooth_cutoff = .false.
      double_count = .false.
      call system_timer('com_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("com_dimer_calc: descriptor object not initialised", error)
      endif

      if (.not. has_property(at, 'mass')) then
         RAISE_ERROR('com_dimer_calc: Atoms has no mass property', error)
      end if
      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      monomer_one_size =size(this%signature_one)
      monomer_two_size =size(this%signature_two)
      dimer_size = monomer_one_size + monomer_two_size

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='com_dimer_calc args_str')) then
            RAISE_ERROR("com_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("com_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      allocate(atomic_index(dimer_size))
      allocate(atomic_index_one(monomer_one_size))
      allocate(atomic_index_two(monomer_two_size))
      allocate(diffs_one(3,monomer_one_size))
      allocate(diffs_two(3,monomer_two_size))
      allocate(weight_one(monomer_one_size))
      allocate(weight_two(monomer_two_size))

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,this%signature_one,associated_to_monomer,this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%monomers_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
      else
        call find_general_monomer(at,monomer_two_index,this%signature_two,associated_to_monomer,this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         call print("WARNING: com_dimer_calc: not all atoms assigned to a monomer, if you have molecules present other than the following, this is OK")
         call print("signature of molecule 1 ")
         call print(this%signature_one)
         call print("signature of molecule 2 ")
         call print(this%signature_two)
      endif

      n_monomer_one = size(monomer_one_index,2)
      n_monomer_two = size(monomer_two_index,2)
      if (n_monomer_one < 1 .or. n_monomer_two < 1) then
        if ( this%strict ) then
          RAISE_ERROR("com_dimer_calc failed to find at least one of the monomer types, try increasing monomer cutoffs", error)
        else
          call print("WARNING: com_dimer_calc failed to find at least one of the monomer types, try increasing monomer cutoffs")
        end if
      end if

      call system_timer('com_dimer_calc: find_monomer_pairs')
      if (this%mpifind) then
         call print("Using find_monomer_pairs_MPI", PRINT_NERD)
         if(associated(atom_mask_pointer)) then
            call find_monomer_pairs_MPI(at,monomer_pairs,com_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,double_count,this%cutoff,error=error,use_com=.true.,atom_mask=atom_mask_pointer)
         else
            call find_monomer_pairs_MPI(at,monomer_pairs,com_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,double_count,this%cutoff,error=error,use_com=.true.)
         end if
      else
         call find_monomer_pairs    (at,monomer_pairs,com_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,double_count,this%cutoff,use_com=.true.,error=error)
      end if
      call system_timer('com_dimer_calc: find_monomer_pairs')

      if ( size(pairs_diffs_map) < 1) then
        if ( this%strict ) then
          RAISE_ERROR("com_dimer_calc did not find any monomer pairs to make a dimer", error)
        else
          call print("WARNING: com_dimer_calc did not find any monomer pairs to make a dimer")
        end if
      end if

      n_descriptors = size(pairs_diffs_map)
      call print("ready to construct "//n_descriptors //" descriptors",PRINT_NERD)
      allocate(descriptor_out%x(n_descriptors))
      loop_descriptor_init: do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(1))
            allocate(descriptor_out%x(i)%ci(dimer_size))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(1,3,dimer_size))
            allocate(descriptor_out%x(i)%ii(dimer_size))
            allocate(descriptor_out%x(i)%pos(3,dimer_size))
            allocate(descriptor_out%x(i)%has_grad_data(dimer_size))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,dimer_size))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo loop_descriptor_init

      if (n_descriptors > 0) then ! only loop over monomers if we actually found any dimers
         i_desc = 0
         loop_monomer_one: do i = 1, n_monomer_one
            if (.not. any(monomer_pairs(1,:) .eq. i)) cycle

            !get indices of monomer and calc internal distances
            atomic_index_one = monomer_one_index(:,i)

            if(associated(atom_mask_pointer)) then
               if(.not. any(atom_mask_pointer(atomic_index_one))) then
                  cycle
               else
                  if(.not. all(atom_mask_pointer(atomic_index_one))) then
                     RAISE_ERROR("com_dimer_calc: atom mask has to encompass either all or none of the atoms of monomer one",error)
                  endif
               endif
            endif

            do i_atomic=1,monomer_one_size
               weight_one(i_atomic) = at%mass(atomic_index_one(i_atomic))
               call print("weight_one("//i_atomic//") = "//weight_one(i_atomic), PRINT_NERD)
            end do
            call print("weight_one = "//weight_one, PRINT_NERD)
            call print("sum(weight_one) = "//sum(weight_one), PRINT_NERD)
            weight_one = weight_one / sum(weight_one)
            call print("weight_one = "//weight_one, PRINT_NERD)

            com_pos_one = centre_of_mass(at, index_list=atomic_index_one)
            !calc atomic positions and shifts relative to mean pos for monomer one, and also distances wrt atom 1
            do i_atomic=1,monomer_one_size
               diffs_one(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_one(i_atomic)),com_pos_one)
            end do

            ! Loop through monomers paired with this one to make dimers
            loop_monomer_pairs: do
               unit_array = maxloc(monomer_pairs(2,:), monomer_pairs(1,:) .eq. i) ! find a monomer paired with i
               this_pair = unit_array(1)

               if (this_pair == 0) exit

               !get indices of monomer two
               j = monomer_pairs(2,this_pair)
               monomer_pairs(:,this_pair) = 0 ! make sure this pair isn't found again
               atomic_index_two = monomer_two_index(:,j)
               atomic_index=(/atomic_index_one,atomic_index_two/)

               do i_atomic=1,monomer_two_size
                  weight_two(i_atomic) = at%mass(atomic_index_two(i_atomic))
               end do
               weight_two = weight_two / sum(weight_two)
               call print("weight_two="//weight_two, PRINT_NERD)

               com_pos_two = centre_of_mass(at, index_list=atomic_index_two)

               ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
               do j_atomic=1,monomer_two_size
                  diffs_two(:,j_atomic) = diff_min_image(at,at%pos(:,atomic_index_two(j_atomic)),com_pos_two)
               end do

               loop_different_shifts: do
                  unit_array = maxloc(pairs_diffs_map, pairs_diffs_map .eq. this_pair) ! find repeats of this pair with different shifts
                  diff_loc = unit_array(1)
                  if (diff_loc == 0) exit

                  i_desc = i_desc + 1
                  
                  diff_one_two = com_pos_diffs(:,diff_loc) ! shift between mean positions of these two monomers
                  pairs_diffs_map(diff_loc) = 0 ! make sure this shifted pair isn't found again

                  ! calculate distance
                  dist = norm(diff_one_two)
                  call print("COM distance = "//dist, PRINT_NERD)

                  calc_descriptor: if(my_do_descriptor) then
                     descriptor_out%x(i_desc)%has_data = .true.
                     if(this%transfer_parameters%do_transfer) then
                        descriptor_out%x(i_desc)%data = transferfunction(dist, this%transfer_parameters)
                     else
                        descriptor_out%x(i_desc)%data = dist
                     end if
                     descriptor_out%x(i_desc)%ci(:) = atomic_index
                     descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp
                     primitive_cutoff = coordination_function(dist,this%cutoff,this%cutoff_transition_width)
                     descriptor_out%x(i_desc)%covariance_cutoff = primitive_cutoff
                  end if calc_descriptor

                  calc_grad_descriptor: if(my_do_grad_descriptor) then !calc grads and update

                     descriptor_out%x(i_desc)%ii(:) = atomic_index
                     descriptor_out%x(i_desc)%has_grad_data(:) = .true.

                     primitive_cutoff_grad = dcoordination_function(dist,this%cutoff,this%cutoff_transition_width)

                     if(this%transfer_parameters%do_transfer) then
                        transdirvec = transferfunction_grad(dist, this%transfer_parameters) * diff_one_two / dist
                     else
                        transdirvec = diff_one_two / dist
                     endif

                     do i_atomic=1,monomer_one_size
                        descriptor_out%x(i_desc)%pos(:,i_atomic) = com_pos_one - diffs_one(:,i_atomic)
                        descriptor_out%x(i_desc)%grad_data(1,:,i_atomic) = - weight_one(i_atomic) * transdirvec ! descriptor wrt atom i_atomic
                        descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = primitive_cutoff_grad * descriptor_out%x(i_desc)%grad_data(1,:,i_atomic)
                     end do

                     do j_atomic=1,monomer_two_size
                        descriptor_out%x(i_desc)%pos(:,monomer_one_size+j_atomic) = com_pos_one + diff_one_two - diffs_two(:,j_atomic)
                        descriptor_out%x(i_desc)%grad_data(1,:,monomer_one_size+j_atomic) = weight_two(j_atomic) * transdirvec ! descriptor wrt atom j_atomic
                        descriptor_out%x(i_desc)%grad_covariance_cutoff(:,monomer_one_size+j_atomic) = primitive_cutoff_grad * descriptor_out%x(i_desc)%grad_data(1,:,monomer_one_size+j_atomic)
                     end do

                  endif calc_grad_descriptor
               enddo loop_different_shifts
            enddo loop_monomer_pairs
         enddo loop_monomer_one
      endif ! (n_descriptors > 0) ... still need to deallocate if no dimers were found:


      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(monomer_pairs)
      deallocate(com_pos_diffs)
      deallocate(pairs_diffs_map)

      deallocate(atomic_index)
      deallocate(atomic_index_one)
      deallocate(atomic_index_two)
      deallocate(weight_one)
      deallocate(weight_two)
      deallocate(diffs_one)
      deallocate(diffs_two)

      deallocate(associated_to_monomer)

      call system_timer('com_dimer_calc')

   endsubroutine com_dimer_calc

   subroutine general_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)

      type(general_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error
      real(dp) :: r_one_two

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      real(dp) :: com_dist, cutoff_grad
      logical :: my_do_descriptor, my_do_grad_descriptor, monomers_identical, use_smooth_cutoff, compound_cutoff
      integer :: d, n_descriptors, n_cross, dimer_size, i, j, k, n, m, i_atomic, j_atomic, start, finish, i_desc, cutoff_pos
      integer :: monomer_one_size, monomer_two_size, n_monomer_one, n_monomer_two, n_products, this_pair, this_shift,diff_loc
      integer, dimension(1) :: unit_array
      real(dp), dimension(3) :: diff_one_two, temp_diff, mean_pos_one, mean_pos_two, transdirvec
      real(dp), dimension(:), allocatable :: dist_vec, primitive_cutoffs, primitive_cutoff_grads
      real(dp), dimension(:,:), allocatable :: interatomic_distances, diffs_one, diffs_two, mean_pos_diffs
      real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
      integer, dimension(3) :: temp_shift, shift_one_two
      real(dp), dimension(:), allocatable :: weight_one, weight_two
      integer, dimension(:), allocatable :: atomic_index, atomic_index_one, atomic_index_two, pairs_diffs_map
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index,  monomer_pairs
      logical, dimension(:), allocatable :: associated_to_monomer


      INIT_ERROR(error)
      use_smooth_cutoff = .false.
      call system_timer('general_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("general_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      monomer_one_size =size(this%signature_one)
      monomer_two_size =size(this%signature_two)
      dimer_size = monomer_one_size + monomer_two_size

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='general_dimer_calc args_str')) then
            RAISE_ERROR("general_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("general_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      if (this%use_com .and. .not. has_property(at, 'mass')) then
         RAISE_ERROR('general_dimer_calc: Atoms has no mass property', error)
      end if

      d = general_dimer_dimensions(this,error)

      compound_cutoff=.True.
      if (count(this%cutoff_contributor) == 1) then
        compound_cutoff=.False.
      else if (count(this%cutoff_contributor) == 0) then
        RAISE_ERROR("general_dimer_calc, initialisation of general dimer did not find a pair of heavy atoms to use for calculating cutoff", error)
      end if

      allocate(dist_vec(d))
      allocate(primitive_cutoffs(d))
      allocate(primitive_cutoff_grads(d))

      allocate(atomic_index(dimer_size))
      allocate(atomic_index_one(monomer_one_size))
      allocate(atomic_index_two(monomer_two_size))
      allocate(diffs_one(3,monomer_one_size))
      allocate(diffs_two(3,monomer_two_size))
      allocate(weight_one(monomer_one_size))
      allocate(weight_two(monomer_two_size))

      allocate(interatomic_vectors(dimer_size,dimer_size,3))
      allocate(interatomic_distances(dimer_size,dimer_size))
      allocate(associated_to_monomer(at%N))
      interatomic_vectors = 0.0_dp
      interatomic_distances = 0.0_dp
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,this%signature_one,associated_to_monomer,this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%monomers_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
      else
        call find_general_monomer(at,monomer_two_index,this%signature_two,associated_to_monomer,this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         call print("WARNING: general_dimer_calc: not all atoms assigned to a monomer, if you have molecules present other than the following, this is OK")
         call print("signature of molecule 1 ")
         call print(this%signature_one)
         call print("signature of molecule 2 ")
         call print(this%signature_two)
      endif

      n_monomer_one = size(monomer_one_index,2)
      n_monomer_two = size(monomer_two_index,2)
      if (n_monomer_one < 1 .or. n_monomer_two < 1) then
        if ( this%strict ) then
          RAISE_ERROR("general_dimer_calc,failed to find at least one of the monomer types, try increasing monomer cutoffs", error)
        else
          call print("WARNING: general_dimer_calc failed to find at least one of the monomer types, try increasing monomer cutoffs")
        end if
      end if

      call system_timer('general_dimer_calc: find_monomer_pairs')
      if (this%mpifind) then
         call print("Using find_monomer_pairs_MPI", PRINT_NERD)
         if(associated(atom_mask_pointer)) then
            call find_monomer_pairs_MPI(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,this%double_count,this%cutoff,error=error,use_com=this%use_com,atom_mask=atom_mask_pointer)
         else
            call find_monomer_pairs_MPI(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,this%double_count,this%cutoff,error=error,use_com=this%use_com)
         end if
      else
         call find_monomer_pairs(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,this%double_count,this%cutoff,error=error,use_com=this%use_com)
      end if
      call system_timer('general_dimer_calc: find_monomer_pairs')

      if ( size(pairs_diffs_map) < 1) then
        if ( this%strict ) then
          RAISE_ERROR("general_dimer_calc did not find any monomer pairs to make a dimer", error)
        else
          call print("WARNING: general_dimer_calc did not find any monomer pairs to make a dimer")
        end if
      end if

      n_descriptors = size(pairs_diffs_map)
      call print("ready to construct "//n_descriptors //" descriptors",PRINT_NERD)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if (my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            allocate(descriptor_out%x(i)%ci(dimer_size))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         end if
         if (my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,dimer_size))
            allocate(descriptor_out%x(i)%ii(dimer_size))
            allocate(descriptor_out%x(i)%pos(3,dimer_size))
            allocate(descriptor_out%x(i)%has_grad_data(dimer_size))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,dimer_size))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         end if
      end do

      has_descriptors: if (n_descriptors > 0) then ! only loop over monomers if we actually found any dimers
         i_desc = 0
         loop_monomer_one: do i = 1, n_monomer_one
            if (.not. any(monomer_pairs(1,:) .eq. i)) cycle

            !get indices of monomer and calc internal distances
            atomic_index_one = monomer_one_index(:,i)

            if (associated(atom_mask_pointer)) then
               if (.not. any(atom_mask_pointer(atomic_index_one))) then
                  cycle
               else
                  if (.not. all(atom_mask_pointer(atomic_index_one))) then
                     RAISE_ERROR("general_dimer_calc: atom mask has to encompass either all or none of the atoms of monomer one",error)
                  end if
               end if
            end if

            if (this%use_com) then
               do i_atomic=1,monomer_one_size
                  weight_one(i_atomic) = at%mass(atomic_index_one(i_atomic))
               end do
               weight_one = weight_one / sum(weight_one)
               mean_pos_one = centre_of_mass(at, index_list=atomic_index_one)
            else
               weight_one = 1.0_dp / monomer_one_size
               mean_pos_one = calc_mean_pos(at,atomic_index_one)
            end if


            !calc atomic positions and shifts relative to mean pos for monomer one, and also distances wrt atom 1
            do i_atomic=1,monomer_one_size
               diffs_one(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_one(i_atomic)),mean_pos_one)
               interatomic_vectors(1,i_atomic,:) = diff_min_image(at,atomic_index_one(1),atomic_index_one(i_atomic))
            end do

            !find other relative positions through vector addition
            do j_atomic=2,monomer_one_size
               do i_atomic=2,j_atomic-1
                  interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:)
               end do
            end do

            !And convert vectors to scalar distances
            do i_atomic=1,monomer_one_size
               do j_atomic=i_atomic+1,monomer_one_size
                  interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
               end do
            end do

            ! Loop through monomers paired with this one to make dimers
            loop_monomer_pairs: do
               unit_array = maxloc(monomer_pairs(2,:), monomer_pairs(1,:) .eq. i) ! find a monomer paired with i
               this_pair = unit_array(1)

               if (this_pair == 0) exit

               !get indices of monomer two
               j = monomer_pairs(2,this_pair)
               monomer_pairs(:,this_pair) = 0 ! make sure this pair isn't found again
               atomic_index_two = monomer_two_index(:,j)
               atomic_index=(/atomic_index_one,atomic_index_two/)

               if (this%use_com) then
                  do j_atomic=1,monomer_two_size
                     weight_two(j_atomic) = at%mass(atomic_index_two(j_atomic))
                  end do
                  weight_two = weight_two / sum(weight_two)
                  mean_pos_two = centre_of_mass(at, index_list=atomic_index_two)
               else
                  weight_two = 1.0_dp / monomer_two_size
                  mean_pos_two = calc_mean_pos(at,atomic_index_two)
               end if

               ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
               do i_atomic=1,monomer_two_size
                  diffs_two(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_two(i_atomic)),mean_pos_two)
                  interatomic_vectors(monomer_one_size+1,monomer_one_size+i_atomic,:) = diff_min_image(at,atomic_index_two(1),atomic_index_two(i_atomic))
               end do

               !find other relative positions through vector addition
               do j_atomic=monomer_one_size+2,dimer_size
                  do i_atomic=monomer_one_size+2,j_atomic-1
                     interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(monomer_one_size+1,j_atomic,:) - interatomic_vectors(monomer_one_size+1,i_atomic,:)
                  end do
               end do

               !And convert vectors to scalar distances
               do i_atomic=monomer_one_size+1,dimer_size
                  do j_atomic=i_atomic+1,dimer_size
                     interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                  end do
               end do

               loop_pair_shifts: do
                  unit_array = maxloc(pairs_diffs_map, pairs_diffs_map .eq. this_pair) ! find repeats of this pair with different shifts
                  diff_loc = unit_array(1)
                  if (diff_loc == 0) exit

                  i_desc = i_desc + 1
                  
                  diff_one_two = mean_pos_diffs(:,diff_loc) ! shift between mean positions of these two monomers
                  pairs_diffs_map(diff_loc) = 0 ! make sure this shifted pair isn't found again

                  com_dist = norm(diff_one_two)

                  ! calculate intermolecular distances, also by vector addition
                  do i_atomic=1,monomer_one_size
                     do j_atomic = monomer_one_size+1,dimer_size
                        interatomic_vectors(i_atomic,j_atomic,:) = diffs_one(:,i_atomic) + diff_one_two - diffs_two(:,j_atomic-monomer_one_size)
                        interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                     end do
                  end do

                  !Now take the whole matrix of scalar distances and combine into 1D array
                  start = 1
                  finish = dimer_size-1
                  do i_atomic=1,dimer_size-1
                     dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:dimer_size)
                     start = finish+1
                     finish=finish + dimer_size-i_atomic-1
                  end do
                  call print("dist vec "//dist_vec,PRINT_NERD)

                  primitive_cutoffs=1.0_dp
                  do k=1,d
                     if (this%cutoff_contributor(k)) then
                        primitive_cutoffs(k) = coordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                     end if
                  end do

                  calc_descriptor: if (my_do_descriptor) then
                     descriptor_out%x(i_desc)%has_data = .true.

                     do_transfer: if (this%transfer_parameters%do_transfer) then
                        do k=1,d
                           if (this%is_intermolecular(k)) then
                              descriptor_out%x(i_desc)%data(k) = transferfunction(dist_vec(k), this%transfer_parameters)
                           else
                              ! don't apply transfer function to intra-molecular (bond) distances
                              descriptor_out%x(i_desc)%data(k) = dist_vec(k)
                           end if
                        end do
                     else
                        descriptor_out%x(i_desc)%data = dist_vec
                     end if do_transfer

                     descriptor_out%x(i_desc)%ci(:) = atomic_index
                     descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp

                     is_com_cutoff: if (this%mpifind) then
                        descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(com_dist, this%cutoff, this%cutoff_transition_width)
                     else

                        is_compound_cutoff: if (compound_cutoff) then
                        ! Covariance cutoff is sum of pairwise products of primitive cutoffs for all *inter*molecular distances excluding H atoms

                           n_products=0
                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                 do m=k+1,d
                                    if (this%cutoff_contributor(m)) then
                                      n_products = n_products + 1
                                      descriptor_out%x(i_desc)%covariance_cutoff = descriptor_out%x(i_desc)%covariance_cutoff + primitive_cutoffs(k)*primitive_cutoffs(m)
                                    end if
                                 end do
                              end if
                           end do
                           ! normalise
                           descriptor_out%x(i_desc)%covariance_cutoff = descriptor_out%x(i_desc)%covariance_cutoff / n_products

                        else ! Covariance cutoff is primitive cutoff, i.e. there is only one pair of heavy atoms

                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                 cutoff_pos=k
                                 descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                                 exit
                              end if
                           end do

                        end if is_compound_cutoff

                     end if is_com_cutoff
                  end if calc_descriptor

                  calc_grad_descriptor: if (my_do_grad_descriptor) then !calc grads and update

                     descriptor_out%x(i_desc)%ii(:) = atomic_index

                     do i_atomic=1,monomer_one_size
                       descriptor_out%x(i_desc)%pos(:,i_atomic) = mean_pos_one - diffs_one(:,i_atomic)
                     end do
                     do i_atomic=1,monomer_two_size
                       descriptor_out%x(i_desc)%pos(:,monomer_one_size+i_atomic) = mean_pos_one + diff_one_two - diffs_two(:,i_atomic)
                     end do

                     call print("DIMER  "//dimer_size,PRINT_NERD)
                     call print('DIMER  Lattice="10.0000000       0.00000000       0.00000000       0.00000000      10.0000000       0.00000000       0.00000000       0.00000000      10.00000" Properties=Z:I:1:pos:R:3',PRINT_NERD)
                     do i_atomic=1,dimer_size
                       call print("DIMER  "//at%Z(atomic_index(i_atomic))//"    "//descriptor_out%x(i_desc)%pos(:,i_atomic),PRINT_NERD)
                     end do


                     !build the grad_data matrix
                     descriptor_out%x(i_desc)%has_grad_data(:) = .true.

                     do k=1,d
                        !get pair of atoms contributing to this component
                        i_atomic = this%component_atoms(k,1)
                        j_atomic = this%component_atoms(k,2)
                        do_transfer_grad: if (this%transfer_parameters%do_transfer .and. this%is_intermolecular(k)) then
                           descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = -transferfunction_grad(dist_vec(k), this%transfer_parameters) * interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic, using a transfer function
                        else
                           descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = -interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic
                        end if do_transfer_grad
                        descriptor_out%x(i_desc)%grad_data(k,:,j_atomic) = -descriptor_out%x(i_desc)%grad_data(k,:,i_atomic)        ! descriptor wrt j_atomic
                     end do

                     is_com_cutoff_grad: if (this%mpifind) then
                        transdirvec = diff_one_two / com_dist

                        cutoff_grad = dcoordination_function(com_dist, this%cutoff, this%cutoff_transition_width)
                        do i_atomic=1,monomer_one_size
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = - weight_one(i_atomic) * transdirvec * cutoff_grad
                        end do
                        do j_atomic=1,monomer_two_size
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,monomer_one_size+j_atomic) = weight_two(j_atomic) * transdirvec * cutoff_grad
                        end do
                     else

                        primitive_cutoff_grads=0.0_dp
                        do k=1,d
                           if (this%cutoff_contributor(k)) then
                             primitive_cutoff_grads(k) = dcoordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                           end if
                        end do

                        is_compound_cutoff_grad: if (compound_cutoff) then

                           do i_atomic=1,dimer_size ! for each atom in the dimer...
                              do k=1,d ! ...iterate over all distance...
                                 if (this%cutoff_contributor(k)) then
                                    do m=k+1,d ! ...pairs involved...
                                       if (this%cutoff_contributor(m)) then
                                          if (any(this%component_atoms(k,:) == i_atomic) .or. any(this%component_atoms(m,:) == i_atomic)) then
                                             descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) &
                                             = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) &
                                                + primitive_cutoffs(m)*primitive_cutoff_grads(k)*descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) &
                                                + primitive_cutoffs(k)*primitive_cutoff_grads(m)*descriptor_out%x(i_desc)%grad_data(m,:,i_atomic)
                                          end if
                                       end if
                                    end do
                                 end if
                              end do
                           end do
                           !normalisation factor
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) / n_products

                        else

                           i_atomic = this%component_atoms(cutoff_pos,1)
                           j_atomic = this%component_atoms(cutoff_pos,2)
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = primitive_cutoff_grads(cutoff_pos) *descriptor_out%x(i_desc)%grad_data(cutoff_pos,:,i_atomic)
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,j_atomic) = primitive_cutoff_grads(cutoff_pos) *descriptor_out%x(i_desc)%grad_data(cutoff_pos,:,j_atomic)

                        end if is_compound_cutoff_grad

                     end if is_com_cutoff_grad
                  end if calc_grad_descriptor

               end do loop_pair_shifts
            end do loop_monomer_pairs
         end do loop_monomer_one
      end if has_descriptors ! (n_descriptors > 0) ... still need to deallocate if no dimers were found:


      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(monomer_pairs)
      deallocate(mean_pos_diffs)
      deallocate(pairs_diffs_map)

      deallocate(dist_vec)
      deallocate(primitive_cutoffs)
      deallocate(primitive_cutoff_grads)

      deallocate(atomic_index)
      deallocate(atomic_index_one)
      deallocate(atomic_index_two)
      deallocate(weight_one)
      deallocate(weight_two)
      deallocate(diffs_one)
      deallocate(diffs_two)

      deallocate(interatomic_vectors)
      deallocate(interatomic_distances)
      deallocate(associated_to_monomer)

      call system_timer('general_dimer_calc')

   end subroutine general_dimer_calc

   subroutine general_trimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)

     type(general_trimer), intent(in) :: this
     type(atoms), intent(in) :: at
     type(descriptor_data), intent(out) :: descriptor_out
     logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
     character(len=*), intent(in), optional :: args_str
     integer, optional, intent(out) :: error
     real(dp) :: r_one_two

     type(Dictionary) :: params
     character(STRING_LENGTH) :: atom_mask_name
     logical :: has_atom_mask_name
     logical, dimension(:), pointer :: atom_mask_pointer

     logical :: my_do_descriptor, my_do_grad_descriptor, one_two_identical, one_three_identical, two_three_identical, &
          use_smooth_cutoff, done_this_monomer, done_this_dimer

     integer :: d, n_descriptors, n_cross, dimer_size, trimer_size, i, j, k, n, m, i_atomic, j_atomic, start, finish, i_desc,this_diff,triplet_pos
     integer :: monomer_one_size, monomer_two_size, monomer_three_size, n_monomer_one, n_monomer_two, n_monomer_three, n_products, contributor_loc
     integer, dimension(1) :: unit_array
     real(dp) :: temp_dist
     real(dp), dimension(3) :: diff_one_two, diff_one_three,diff_two_three,mean_pos_one,mean_pos_two,mean_pos_three
     real(dp), dimension(3) :: grad_cut_one, grad_cut_two, grad_cut_three
     real(dp) :: dist_one_two, dist_one_three, dist_two_three, cut12, cut13, cut23, dcut12, dcut13, dcut23
     real(dp), dimension(:), allocatable :: dist_vec, primitive_cutoffs, primitive_cutoff_grads
     real(dp), dimension(:,:), allocatable :: interatomic_distances,diffs_one,diffs_two,diffs_three,triplets_diffs
     real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
     integer, dimension(3) :: temp_shift, shift_one_two,shift_one_three
     integer, dimension(:), allocatable :: atomic_index_dimer, atomic_index_trimer, atomic_index_one, atomic_index_two, atomic_index_three,triplets_diffs_map
     integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index, monomer_three_index, monomer_triplets
     logical, dimension(:), allocatable :: associated_to_monomer
     logical :: double_count

     INIT_ERROR(error)
     use_smooth_cutoff = .false.
     double_count=.false.
     call system_timer('general_trimer_calc')

     if(.not. this%initialised) then
        RAISE_ERROR("general_trimer_calc: descriptor object not initialised", error)
     endif

     my_do_descriptor = optional_default(.false., do_descriptor)
     my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

     monomer_one_size =size(this%signature_one)
     monomer_two_size =size(this%signature_two)
     monomer_three_size =size(this%signature_three)
     dimer_size = monomer_one_size + monomer_two_size
     trimer_size = monomer_one_size + monomer_two_size + monomer_three_size

     if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

     call finalise(descriptor_out)

     atom_mask_pointer => null()
     if(present(args_str)) then
        call initialise(params)
        
        call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
        help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
        "calculated.")

        if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='general_trimer_calc args_str')) then
           RAISE_ERROR("general_trimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
        endif
        
        call finalise(params)

        if( has_atom_mask_name ) then
           if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
              RAISE_ERROR("general_trimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
           endif
        else
           atom_mask_pointer => null()
        endif

     endif

     d = general_trimer_dimensions(this,error)

     allocate(dist_vec(d))
     allocate(primitive_cutoffs(d))
     allocate(primitive_cutoff_grads(d))

     allocate(atomic_index_one(monomer_one_size))
     allocate(atomic_index_two(monomer_two_size))
     allocate(atomic_index_three(monomer_three_size))
     allocate(atomic_index_dimer(dimer_size))
     allocate(atomic_index_trimer(trimer_size))

     allocate(interatomic_vectors(trimer_size,trimer_size,3))
     allocate(interatomic_distances(trimer_size,trimer_size))
     allocate(associated_to_monomer(at%N))
     allocate(diffs_one(3,monomer_one_size))
     allocate(diffs_two(3,monomer_two_size))
     allocate(diffs_three(3,monomer_three_size))

     interatomic_vectors = 0.0_dp
     interatomic_distances = 0.0_dp
     associated_to_monomer=.false.

     call find_general_monomer(at,monomer_one_index,this%signature_one,associated_to_monomer,this%monomer_one_cutoff,this%atom_ordercheck,error)
     if (this%one_two_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
     else
        call find_general_monomer(at,monomer_two_index,this%signature_two,associated_to_monomer,this%monomer_two_cutoff,this%atom_ordercheck,error)
     end if
     if (this%one_three_identical) then
        allocate(monomer_three_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_three_index = monomer_one_index
     else if (this%two_three_identical) then
        allocate(monomer_three_index(size(monomer_two_index,1),size(monomer_two_index,2)))
        monomer_three_index = monomer_two_index
     else
        call find_general_monomer(at,monomer_three_index,this%signature_three,associated_to_monomer,this%monomer_three_cutoff,this%atom_ordercheck,error)
     end if

     if(.not. all(associated_to_monomer)) then
        if(this%strict) then
           RAISE_ERROR("general_trimer_calc: not all atoms assigned to a monomer", error)
        else
           call print("WARNING: general_trimer_calc: not all atoms assigned to a monomer")
        endif
     endif

     n_monomer_one = size(monomer_one_index,2)
     n_monomer_two = size(monomer_two_index,2)
     n_monomer_three = size(monomer_three_index,2)

     if (this%use_com) then
        RAISE_ERROR("general_trimer_calc: use_com=T not implemented yet", error)
     end if
     call system_timer('general_trimer_calc: find_monomer_triplets')
     if(this%mpifind) then
        call print("Using find_monomer_triplets_MPI", PRINT_NERD)
        if(associated(atom_mask_pointer)) then
           call find_monomer_triplets_MPI(at,monomer_triplets,triplets_diffs,triplets_diffs_map,monomer_one_index,monomer_two_index,monomer_three_index,this%one_two_identical,this%one_three_identical,this%two_three_identical,this%cutoff,error,use_com=.false.,atom_mask=atom_mask_pointer)
        else
           call find_monomer_triplets_MPI(at,monomer_triplets,triplets_diffs,triplets_diffs_map,monomer_one_index,monomer_two_index,monomer_three_index,this%one_two_identical,this%one_three_identical,this%two_three_identical,this%cutoff,error,use_com=.false.)
        end if
     else
       call find_monomer_triplets(at,monomer_triplets,triplets_diffs,triplets_diffs_map,monomer_one_index,monomer_two_index,monomer_three_index,this%one_two_identical,this%one_three_identical,this%two_three_identical,this%cutoff,error)
     end if
     call system_timer('general_trimer_calc: find_monomer_triplets')
     call print("monomer_triplets("//size(monomer_triplets,1)//", "//size(monomer_triplets,2)//")", PRINT_NERD)
     call print("triplets_diffs("//size(triplets_diffs,1)//", "//size(triplets_diffs,2)//")", PRINT_NERD)
     call print("triplets_diffs_map("//size(triplets_diffs_map)//")", PRINT_NERD)

     !call print("monomer_triplets",PRINT_NERD)
     !call print(monomer_triplets,PRINT_NERD)

     n_descriptors = size(triplets_diffs_map)
     call print("ready to construct "//n_descriptors //" descriptors",PRINT_NERD)
     allocate(descriptor_out%x(n_descriptors))
     do i = 1, n_descriptors
        if(my_do_descriptor) then
           allocate(descriptor_out%x(i)%data(d))
           allocate(descriptor_out%x(i)%ci(trimer_size))
           descriptor_out%x(i)%data = 0.0_dp
           descriptor_out%x(i)%ci = 0
           descriptor_out%x(i)%has_data = .false.
           descriptor_out%x(i)%covariance_cutoff = 1.0_dp
        endif
        if(my_do_grad_descriptor) then
           allocate(descriptor_out%x(i)%grad_data(d,3,trimer_size))
           allocate(descriptor_out%x(i)%ii(trimer_size))
           allocate(descriptor_out%x(i)%pos(3,trimer_size))
           allocate(descriptor_out%x(i)%has_grad_data(trimer_size))
           descriptor_out%x(i)%grad_data = 0.0_dp
           descriptor_out%x(i)%ii = 0
           descriptor_out%x(i)%pos = 0.0_dp
           descriptor_out%x(i)%has_grad_data = .false.

           allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,trimer_size))
           descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
        endif
     enddo

      has_descriptors: if (n_descriptors > 0) then ! only loop over monomers if we actually found any trimers
         i_desc = 0
         loop_monomer_one: do i = 1, n_monomer_one
            if (.not. any(monomer_triplets(1,:) .eq. i)) cycle
            done_this_monomer = .false.
            !get indices of monomer and calc internal distances
            atomic_index_one = monomer_one_index(:,i) !store the indices of atoms in this monomer

            if(associated(atom_mask_pointer)) then
               if(.not. any(atom_mask_pointer(atomic_index_one))) then
                  cycle
               else
                  if(.not. all(atom_mask_pointer(atomic_index_one))) then
                     RAISE_ERROR("general_trimer_calc: atom mask has to encompass either all or none of the atoms of monomer one",error)
                  endif
               endif
            endif

            mean_pos_one = calc_mean_pos(at,atomic_index_one)

            !calc atomic positions and shifts relative to mean pos for monomer one, and also distances wrt atom 1
            do i_atomic=1,monomer_one_size
               diffs_one(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_one(i_atomic)),mean_pos_one)
               interatomic_vectors(1,i_atomic,:) = diff_min_image(at,atomic_index_one(1),atomic_index_one(i_atomic))
            end do

            !find other relative positions through vector addition
            do j_atomic=2,monomer_one_size
               do i_atomic=2,j_atomic-1
                  interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:)
               end do
            end do

            !Now convert vectors to scalar distances
            do i_atomic=1,monomer_one_size
               do j_atomic=i_atomic+1,monomer_one_size
                  interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
               end do
            end do

            ! Loop through monomers paired with this one to make dimers
            loop_monomer_pairs: do
               unit_array = maxloc(monomer_triplets(2,:), monomer_triplets(1,:) .eq. i) ! find a monomer paired with i
               if (all(unit_array .eq. 0)) exit

               !get indices of monomer two
               j = monomer_triplets(2,unit_array(1))
               atomic_index_two = monomer_two_index(:,j)
               atomic_index_dimer=(/atomic_index_one,atomic_index_two/)

               mean_pos_two = calc_mean_pos(at,atomic_index_two)

               ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
               do i_atomic=1,monomer_two_size
                  diffs_two(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_two(i_atomic)),mean_pos_two)
                  interatomic_vectors(monomer_one_size+1,monomer_one_size+i_atomic,:) = diff_min_image(at,atomic_index_two(1),atomic_index_two(i_atomic))
               end do

               !find other relative positions through vector addition
               do j_atomic=monomer_one_size+2,dimer_size
                  do i_atomic=monomer_one_size+2,j_atomic-1
                     interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(monomer_one_size+1,j_atomic,:) - interatomic_vectors(monomer_one_size+1,i_atomic,:)
                  end do
               end do

               !And convert vectors to scalar distances
               do i_atomic=monomer_one_size+1,dimer_size
                  do j_atomic=i_atomic+1,dimer_size
                     interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                  end do
               end do

               !! Now make trimers based on this dimer
               loop_triplets: do
                  unit_array = maxloc(monomer_triplets(3,:), monomer_triplets(1,:) .eq. i .and. monomer_triplets(2,:) .eq. j ) ! look for  a triplet
                  if (all(unit_array .eq. 0)) exit
                  triplet_pos=unit_array(1)

                  !get indices of monomer three
                  k = monomer_triplets(3,triplet_pos)
                  monomer_triplets(:,triplet_pos) = 0 ! make sure this triplet isn't found again
                  atomic_index_three = monomer_three_index(:,k)
                  atomic_index_trimer=(/atomic_index_dimer,atomic_index_three/)
                  mean_pos_three = calc_mean_pos(at,atomic_index_three)
                  ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
                  do i_atomic=1,monomer_three_size
                     diffs_three(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_three(i_atomic)),mean_pos_three)
                     interatomic_vectors(dimer_size+1,dimer_size+i_atomic,:) = diff_min_image(at,atomic_index_three(1),atomic_index_three(i_atomic))
                  end do

                  !find other relative positions through vector addition
                  do j_atomic=dimer_size+2,trimer_size
                     do i_atomic=dimer_size+2,j_atomic-1
                        interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(dimer_size+1,j_atomic,:) -interatomic_vectors(dimer_size+1,i_atomic,:)
                     end do
                  end do

                  !Now convert vectors to scalar distances
                  do i_atomic=dimer_size+1,trimer_size
                     do j_atomic=i_atomic+1,trimer_size
                        interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                     end do
                  end do

                  ! Loop over all trimers which can be created from these three monomers, i.e. with different periodic images
                  loop_trimers: do
                     unit_array = maxloc(triplets_diffs_map, triplets_diffs_map .eq. triplet_pos)
                     this_diff = unit_array(1)
                     if (this_diff == 0) exit

                     i_desc = i_desc + 1

                     diff_one_two = triplets_diffs(1:3,this_diff)
                     diff_one_three = triplets_diffs(4:6,this_diff)
                     diff_two_three = diff_one_three - diff_one_two
                     triplets_diffs_map(this_diff)=0 ! make sure this shifted triplet isn't found again
                     dist_one_two = norm(diff_one_two)
                     dist_one_three = norm(diff_one_three)
                     dist_two_three = norm(diff_two_three)

                     ! calculate intermolecular distances, monomers one and two
                     do i_atomic=1,monomer_one_size
                        do j_atomic = monomer_one_size+1,dimer_size
                           interatomic_vectors(i_atomic,j_atomic,:) = diffs_one(:,i_atomic) + diff_one_two - diffs_two(:,j_atomic-monomer_one_size)
                           interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                        end do
                     end do

                     ! calculate intermolecular distances, monomers one and three
                     do i_atomic=1,monomer_one_size
                        do j_atomic = dimer_size+1,trimer_size
                           interatomic_vectors(i_atomic,j_atomic,:) = diffs_one(:,i_atomic) + diff_one_three - diffs_three(:,j_atomic-dimer_size)
                           interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                        end do
                     end do

                     ! calculate intermolecular distances, monomers two and three
                     do i_atomic=monomer_one_size+1,dimer_size
                        do j_atomic = dimer_size+1,trimer_size
                           interatomic_vectors(i_atomic,j_atomic,:) = diffs_two(:,i_atomic-monomer_one_size) + diff_two_three - diffs_three(:,j_atomic-dimer_size)
                           interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                        end do
                     end do

                     !Now take the whole matrix of scalar distances and combine into 1D array
                     start = 1
                     finish=trimer_size-1
                     do i_atomic=1,trimer_size-1
                        dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:trimer_size)
                        start = finish+1
                        finish=finish + trimer_size-i_atomic-1
                     end do
                     call print( "list of distances: "//dist_vec,PRINT_NERD)

                     calc_descriptor: if(my_do_descriptor) then
                        descriptor_out%x(i_desc)%has_data = .true.
                        descriptor_out%x(i_desc)%data = dist_vec
                        descriptor_out%x(i_desc)%ci(:) = atomic_index_trimer

                        cutoff_type: if (this%mpifind) then
                           cut12 = coordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
                           cut13 = coordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
                           cut23 = coordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
                           descriptor_out%x(i_desc)%covariance_cutoff = (cut12*cut13 + cut12*cut23 + cut13*cut23)/3.0_dp
                        else
                           descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp

                           primitive_cutoffs=1.0_dp

                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                 primitive_cutoffs(k) = coordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                              end if
                           end do
                           n_products=0
                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                 do m=k+1,d
                                    if (this%cutoff_contributor(m)) then
                                       n_products=n_products+1
                                       descriptor_out%x(i_desc)%covariance_cutoff =  descriptor_out%x(i_desc)%covariance_cutoff + primitive_cutoffs(k)*primitive_cutoffs(m)
                                    end if
                                 end do
                              end if
                           end do
                           ! normalise
                           descriptor_out%x(i_desc)%covariance_cutoff =  descriptor_out%x(i_desc)%covariance_cutoff / n_products
                        end if cutoff_type
                     end if calc_descriptor

                     calc_grad_descriptor: if(my_do_grad_descriptor) then !calc grads and update

                        descriptor_out%x(i_desc)%ii(:) = atomic_index_trimer
                        do i_atomic=1,monomer_one_size
                           descriptor_out%x(i_desc)%pos(:,i_atomic) = mean_pos_one - diffs_one(:,i_atomic)
                        end do
                        do i_atomic=1,monomer_two_size
                           descriptor_out%x(i_desc)%pos(:,monomer_one_size+i_atomic) = mean_pos_one + diff_one_two - diffs_two(:,i_atomic)
                        end do
                        do i_atomic=1,monomer_three_size
                           descriptor_out%x(i_desc)%pos(:,dimer_size+i_atomic) = mean_pos_one + diff_one_three - diffs_three(:,i_atomic)
                        end do

                        call print("TRIM  9",PRINT_NERD)
                        call print('TRIM  Lattice="10.0000000       0.00000000       0.00000000       0.00000000      10.0000000       0.00000000       0.00000000       0.00000000      10.00000" Properties=Z:I:1:pos:R:3',PRINT_NERD)
                        do i_atomic=1,trimer_size
                        call print("TRIM  "//at%Z(atomic_index_trimer(i_atomic))//"    "//descriptor_out%x(i_desc)%pos(:,i_atomic),PRINT_NERD)
                        end do

                        !build the grad_data matrix
                        descriptor_out%x(i_desc)%has_grad_data(:) = .true.
                        do k=1,d
                           i_atomic = this%component_atoms(k,1)
                           j_atomic = this%component_atoms(k,2)
                           descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = -interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic
                           descriptor_out%x(i_desc)%grad_data(k,:,j_atomic) = -descriptor_out%x(i_desc)%grad_data(k,:,i_atomic)        ! descriptor wrt j_atomic
                        end do

                        cutoff_type_grad: if (this%mpifind) then
                           cut12 = coordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
                           cut13 = coordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
                           cut23 = coordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
                           dcut12 = dcoordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
                           dcut13 = dcoordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
                           dcut23 = dcoordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
                           !descriptor_out%x(i_desc)%covariance_cutoff = (cut12*cut13 + cut12*cut23 + cut13*cut23)/3.0_dp
                           grad_cut_one = (- diff_one_two / dist_one_two * dcut12 * (cut13+cut23) &
                                           - diff_one_three / dist_one_three * dcut13 * (cut12+cut23)) / 3.0_dp
                           grad_cut_two = (  diff_one_two / dist_one_two * dcut12 * (cut13+cut23) &
                                           - diff_two_three / dist_two_three * dcut23 * (cut12+cut13)) / 3.0_dp
                           grad_cut_three = (diff_one_three / dist_one_three * dcut13 * (cut12+cut23) &
                                           + diff_two_three / dist_two_three * dcut23 * (cut12+cut13)) / 3.0_dp
                           do i_atomic=1,monomer_one_size
                              descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = grad_cut_one / real(monomer_one_size,dp)
                           end do
                           do i_atomic=1,monomer_two_size
                              descriptor_out%x(i_desc)%grad_covariance_cutoff(:,monomer_one_size+i_atomic) = grad_cut_two / real(monomer_two_size,dp)
                           end do
                           do i_atomic=1,monomer_three_size
                              descriptor_out%x(i_desc)%grad_covariance_cutoff(:,dimer_size+i_atomic) = grad_cut_three / real(monomer_three_size,dp)
                           end do
                        else
                           primitive_cutoff_grads=0.0_dp
                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                primitive_cutoff_grads(k) = dcoordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                              end if
                           end do

                           do i_atomic=1,trimer_size
                             do k=1,d
                                if (this%cutoff_contributor(k) ) then
                                  do m=k+1,d
                                    if (this%cutoff_contributor(m)) then
                                      if (any(this%component_atoms(k,:) == i_atomic) .or. any(this%component_atoms(m,:) == i_atomic)) then
                                      descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) &
                                        + primitive_cutoffs(m)*primitive_cutoff_grads(k)*descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) &
                                        + primitive_cutoffs(k)*primitive_cutoff_grads(m)*descriptor_out%x(i_desc)%grad_data(m,:,i_atomic)
                                      end if
                                    end if
                                  end do
                                end if
                             end do
                           end do

                           !normalisation factor
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) / n_products
                        end if cutoff_type_grad

                     end if calc_grad_descriptor

                  end do loop_trimers
               end do loop_triplets
            end do loop_monomer_pairs
         end do loop_monomer_one
      end if has_descriptors ! (n_descriptors > 0) ... still need to deallocate if no trimers were found:

     deallocate(monomer_one_index)
     deallocate(monomer_two_index)
     deallocate(monomer_three_index)
     if (allocated(monomer_triplets)) deallocate(monomer_triplets)
     if (allocated(triplets_diffs)) deallocate(triplets_diffs)
     if (allocated(triplets_diffs_map)) deallocate(triplets_diffs_map)

     deallocate(dist_vec)
     deallocate(primitive_cutoffs)
     deallocate(primitive_cutoff_grads)

     deallocate(atomic_index_dimer)
     deallocate(atomic_index_trimer)
     deallocate(atomic_index_one)
     deallocate(atomic_index_two)
     deallocate(atomic_index_three)
     deallocate(diffs_one)
     deallocate(diffs_two)
     deallocate(diffs_three)

     deallocate(interatomic_vectors)
     deallocate(interatomic_distances)
     deallocate(associated_to_monomer)


     call system_timer('general_trimer_calc')

   endsubroutine general_trimer_calc

   subroutine rdf_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(rdf), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, i_n, l_n_neighbours, i_desc, n_descriptors, n_cross
      integer, dimension(3) :: shift
      real(dp) :: r_ij, f_cut, df_cut
      real(dp), dimension(3) :: u_ij
      real(dp), dimension(:), allocatable :: rdf_ij

      INIT_ERROR(error)

      call system_timer('rdf_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("rdf_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='rdf_calc args_str')) then
            RAISE_ERROR("rdf_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("rdf_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      call finalise(descriptor_out)

      d = rdf_dimensions(this,error)
      allocate(rdf_ij(d))

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc = i_desc + 1
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%has_data = .false.

            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         i_n = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)

            if( r_ij >= this%cutoff ) cycle
            i_n = i_n + 1

            rdf_ij = exp( -0.5_dp * (r_ij - this%r_gauss)**2 / this%w_gauss**2 )
            f_cut = coordination_function(r_ij,this%cutoff,this%transition_width)

            if(my_do_descriptor) &
               descriptor_out%x(i_desc)%data = descriptor_out%x(i_desc)%data + rdf_ij * f_cut

            if(my_do_grad_descriptor) then
               df_cut = dcoordination_function(r_ij,this%cutoff,this%transition_width)

               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.

               descriptor_out%x(i_desc)%grad_data(:,:,i_n) = ( - ( rdf_ij * (r_ij - this%r_gauss) / this%w_gauss**2 ) * f_cut + rdf_ij * df_cut ) .outer. u_ij
               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,i_n)
            endif
         enddo
      enddo

      if(allocated(rdf_ij)) deallocate(rdf_ij)

      call system_timer('rdf_calc')

   endsubroutine rdf_calc

   subroutine as_distance_2b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(as_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(descriptor) :: my_coordination
      type(descriptor_data) :: descriptor_coordination

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zi1, Zi2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m, n_neighbours_coordination_i, n_neighbours_coordination_ij
      integer, dimension(3) :: shift
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, cos_jik, f_cut_i, f_cut_j, f_cut_ij, f_cut_ik, f_cut_jk, f_cut_as_i, f_cut_as_j, rho_i, rho_j
      real(dp), dimension(3) :: u_ij, u_ik, u_jk

      INIT_ERROR(error)
      call system_timer('as_distance_2b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("as_distance_2b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='as_distance_2b_calc args_str')) then
            RAISE_ERROR("as_distance_2b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("as_distance_2b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("as_distance_2b_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = as_distance_2b_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)


      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance=r_ij, cosines=u_ij)

            if(r_ij > this%max_cutoff .or. r_ij < this%min_cutoff) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            rho_i = 0.0_dp
            f_cut_i = 0.0_dp

            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance=r_ik, cosines=u_ik)

               if(r_ik > this%coordination_cutoff) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               f_cut_ik = coordination_function(r_ik,this%coordination_cutoff,this%coordination_transition_width)

               f_cut_i = f_cut_i + f_cut_ik
               rho_i = rho_i + 0.5_dp * ( erf(cos_ijk/this%overlap_alpha) + 1.0_dp ) * f_cut_ik**2
            enddo

            rho_i = rho_i / f_cut_i

            if(rho_i > this%as_cutoff) cycle

            rho_j = 0.0_dp
            f_cut_j = 0.0_dp

            do m = 1, n_neighbours(at,j)
               k = neighbour(at, j, m, distance=r_jk, cosines=u_jk)

               if(r_jk > this%coordination_cutoff) cycle

               cos_jik = dot_product(-u_ij,u_jk)
               f_cut_jk = coordination_function(r_jk,this%coordination_cutoff,this%coordination_transition_width)

               f_cut_j = f_cut_j + f_cut_jk
               rho_j = rho_j + 0.5_dp * ( erf(cos_jik/this%overlap_alpha) + 1.0_dp ) * f_cut_jk**2
            enddo

            if(rho_j > this%as_cutoff) cycle
            ! all three conditions fulfilled: pair within lower and upper cutoff, asymmetricity lower than threshold

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc)%data(d))
               descriptor_out%x(i_desc)%data = 0.0_dp
               allocate(descriptor_out%x(i_desc)%ci(2))
               descriptor_out%x(i_desc)%has_data = .false.
            endif

            if(my_do_grad_descriptor) then
               n_neighbours_coordination_ij = n_neighbours(at,i,max_dist=this%coordination_cutoff) + &
               n_neighbours(at,j,max_dist=this%coordination_cutoff) + 2

               allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%ii(0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%pos(3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%has_grad_data(0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_data = 0.0_dp
               descriptor_out%x(i_desc)%ii = 0
               descriptor_out%x(i_desc)%pos = 0.0_dp
               descriptor_out%x(i_desc)%has_grad_data = .false.

               allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
            endif
         enddo
      enddo

      i_desc = 0
      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)

            if(r_ij > this%max_cutoff .or. r_ij < this%min_cutoff) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            rho_i = 0.0_dp
            f_cut_i = 0.0_dp

            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance=r_ik, cosines=u_ik)

               if(r_ik > this%coordination_cutoff) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               f_cut_ik = coordination_function(r_ik,this%coordination_cutoff,this%coordination_transition_width)

               f_cut_i = f_cut_i + f_cut_ik
               rho_i = rho_i + 0.5_dp * ( erf(cos_ijk/this%overlap_alpha) + 1.0_dp ) * f_cut_ik**2
            enddo

            rho_i = rho_i / f_cut_i

            if(rho_i > this%as_cutoff) cycle

            rho_j = 0.0_dp
            f_cut_j = 0.0_dp

            do m = 1, n_neighbours(at,j)
               k = neighbour(at, j, m, distance=r_jk, cosines=u_jk)

               if(r_jk > this%coordination_cutoff) cycle

               cos_jik = dot_product(-u_ij,u_jk)
               f_cut_jk = coordination_function(r_jk,this%coordination_cutoff,this%coordination_transition_width)

               f_cut_j = f_cut_j + f_cut_jk
               rho_j = rho_j + 0.5_dp * ( erf(cos_jik/this%overlap_alpha) + 1.0_dp ) * f_cut_jk**2
            enddo

            if(rho_j > this%as_cutoff) cycle
            ! all three conditions fulfilled: pair within lower and upper cutoff, asymmetricity lower than threshold

            i_desc = i_desc + 1

            f_cut_ij = coordination_function(r_ij,this%max_cutoff,this%max_transition_width,this%min_cutoff,this%min_transition_width)
            f_cut_as_i = coordination_function(rho_i,this%as_cutoff, this%as_transition_width)
            f_cut_as_j = coordination_function(rho_j,this%as_cutoff, this%as_transition_width)

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(1:2) = (/i,j/)
               
               descriptor_out%x(i_desc)%has_data = .true.

               descriptor_out%x(i_desc)%data(1) = r_ij
               descriptor_out%x(i_desc)%data(2) = f_cut_i + f_cut_j
               descriptor_out%x(i_desc)%data(3) = (f_cut_i - f_cut_j)**2

               descriptor_out%x(i_desc)%covariance_cutoff = f_cut_ij * f_cut_as_i * f_cut_as_j
            endif
            if(my_do_grad_descriptor) then
               n_neighbours_coordination_i = n_neighbours(at,i,max_dist=this%coordination_cutoff)

               descriptor_out%x(i_desc)%ii(0) = i
               descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
               descriptor_out%x(i_desc)%has_grad_data(0) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,0) = -u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = -dcoordination_function(r_ij,this%coordination_cutoff,this%coordination_transition_width)*u_ij

               !descriptor_out%x(i_desc)%ii(1) = j
               !descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift)
               !descriptor_out%x(i_desc)%has_grad_data(1) = .true.
               !descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij(:)
               !descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0)

               !descriptor_out%x(i_desc)%ii(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%ii(:)
               !descriptor_out%x(i_desc)%pos(:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%pos(:,:)
               !descriptor_out%x(i_desc)%has_grad_data(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%has_grad_data(:)
               !descriptor_out%x(i_desc)%grad_data(2,:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%grad_data(1,:,:)
               !descriptor_out%x(i_desc)%grad_data(3,:,2:n_neighbours_coordination_i+2) = 2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
               !   descriptor_coordination%x(i)%grad_data(1,:,:)

               !descriptor_out%x(i_desc)%ii(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%ii(:)
               !descriptor_out%x(i_desc)%pos(:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%pos(:,:)
               !descriptor_out%x(i_desc)%has_grad_data(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%has_grad_data(:)
               !descriptor_out%x(i_desc)%grad_data(2,:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%grad_data(1,:,:)
               !descriptor_out%x(i_desc)%grad_data(3,:,n_neighbours_coordination_i+3:) = -2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
               !   descriptor_coordination%x(j)%grad_data(1,:,:)

            endif
         enddo
      enddo

      call finalise(my_coordination)
      call finalise(descriptor_coordination)

      call system_timer('as_distance_2b_calc')

   endsubroutine as_distance_2b_calc

   subroutine molecule_lo_d_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(molecule_lo_d), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, i, j, i_atomic, j_atomic, k, start, finish, molecule_size, i_component, i_desc
      integer, dimension(3) :: temp_shift
      real(dp), dimension(:), allocatable :: dist_vec
      real(dp), dimension(:,:), allocatable :: interatomic_distances
      real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
      integer, dimension(:), allocatable :: atomic_index
      integer, dimension(:,:), allocatable :: shifts
      logical, dimension(:), allocatable :: associated_to_molecule


      INIT_ERROR(error)

      call system_timer('molecule_lo_d_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("molecule_lo_d_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='molecule_lo_d_calc args_str')) then
            RAISE_ERROR("molecule_lo_d_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("molecule_lo_d_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("molecule_lo_d_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      molecule_size=this%n_atoms
      d = molecule_lo_d_dimensions(this,error)

      if (this%atom_ordercheck) then
        do i=1,molecule_size
          if (this%template_atoms%Z(i) /= at%Z(i)) then
            call print("atoms in all input frams should be in the same order as this template : ")
            call print(this%template_atoms)
            RAISE_ERROR("molecule_lo_d_calc: atoms not in same order as in template used to teach", error)
          end if
        end do
      end if

      allocate(shifts(molecule_size,3))
      allocate(dist_vec(this%max_dimension))
      allocate(atomic_index(molecule_size))

      allocate(associated_to_molecule(at%N))
      associated_to_molecule=.True.
      do i=1,molecule_size
        atomic_index(i) = i
      end do


      allocate(interatomic_vectors(molecule_size,molecule_size,3))
      allocate(interatomic_distances(molecule_size,molecule_size))
      interatomic_vectors = 0.0_dp
      interatomic_distances = 0.0_dp

!      call find_general_monomer(at,monomer_index,this%signature,associated_to_monomer,this%cutoff,this%atom_ordercheck,error)

      if(.not. all(associated_to_molecule)) then
         RAISE_ERROR("molecule_lo_d_calc: not all atoms assigned to a monomer", error)
      endif

      n_descriptors = 1 ! currently no support for finding multiple nontrivial molecules, otherwise would be size(monomer_index,2)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            allocate(descriptor_out%x(i)%ci(molecule_size))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,molecule_size))
            allocate(descriptor_out%x(i)%ii(molecule_size))
            allocate(descriptor_out%x(i)%pos(3,molecule_size))
            allocate(descriptor_out%x(i)%has_grad_data(molecule_size))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,molecule_size))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo


      do i_desc = 1, n_descriptors

         !atomic_index = monomer_index(:,i) !for fixed ordering don't need this
         ! for now calculate all O(N^2) distances and just pick out the ones we want

         !calc all positions relative to atom 1
         do i_atomic=2,molecule_size
           temp_shift=0
           interatomic_vectors(1,i_atomic,:) = diff_min_image(at,atomic_index(1),atomic_index(i_atomic),shift=temp_shift)
           shifts(i_atomic,:) = temp_shift
         end do

         !find other relative positions through vector addition
         do j_atomic=2,molecule_size
           do i_atomic=2,j_atomic-1
             interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:) 
           end do
         end do

         !Now convert vectors to scalar distances
         do i_atomic=1,molecule_size
           do j_atomic=i_atomic+1,molecule_size
             interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
           end do
         end do       

         !and convert this NxN matrix into the required vector length N(N-1)/2
         start = 1
         finish = molecule_size-1
         do i_atomic=1,molecule_size-1
           dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:molecule_size)  
           start = finish+1
           finish=finish + molecule_size-i_atomic-1
         end do



         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(:) = atomic_index
            descriptor_out%x(i_desc)%has_data = .true.
            if(this%desctype == 0) then
               descriptor_out%x(i_desc)%data = dist_vec(this%included_components) ! here's where we pick out the components
            else if(this%desctype == 1) then
               descriptor_out%x(i_desc)%data = 1.0_dp / dist_vec(this%included_components)
            else if(this%desctype == 2) then
               do k=1,d
                  i_component = this%included_components(k)
                  i_atomic = this%component_atoms(i_component,1)
                  j_atomic = this%component_atoms(i_component,2)
                  descriptor_out%x(i_desc)%data(k) = real(at%Z(i_atomic)*at%Z(j_atomic)) / dist_vec(i_component)
               end do
            else if(this%desctype == 3) then
               descriptor_out%x(i_desc)%data = exp(-dist_vec(this%included_components))
            else if(this%desctype < 0) then
               descriptor_out%x(i_desc)%data = dist_vec(this%included_components)**this%desctype
            else
               RAISE_ERROR("molecule_lo_d_calc: not implemented descriptor type",error)
            end if
         endif
!call print(descriptor_out%x(i_desc)%data)
!!$
!!$        do i_component=1,d
!!$            write(*,*) 'component : ',i_component
!!$            write(*,*) 'value : ',descriptor_out%x(i_desc)%data(i_component)
!!$            write(*,*) 'atoms : ',this%component_atoms(this%included_components(i_component),:)
!!$            write(*,*) 'dist_mat_entry : ',interatomic_distances(this%component_atoms(this%included_components(i_component),1),this%component_atoms(this%included_components(i_component),2))
!!$        end do

         if(my_do_grad_descriptor) then

            descriptor_out%x(i_desc)%ii(:) = atomic_index
            descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,atomic_index(1))
            do i_atomic =2,molecule_size
              descriptor_out%x(i_desc)%pos(:,i_atomic) = at%pos(:,atomic_index(i_atomic)) + matmul(at%lattice,shifts(i_atomic,:))
            end do

            !build the grad_data matrix
            descriptor_out%x(i_desc)%has_grad_data(:) = .true.


            do k=1,d
              !get pair of atoms contributing to this component
              i_component = this%included_components(k)
              i_atomic = this%component_atoms(i_component,1)
              j_atomic = this%component_atoms(i_component,2)

              if(this%desctype == 0) then
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = -interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic 
              else if(this%desctype == 1) then
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)**3  ! descriptor wrt atom i_atomic 
              else if(this%desctype == 2) then
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = real(at%Z(i_atomic)*at%Z(j_atomic)) * interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)**3  ! descriptor wrt atom i_atomic
              else if(this%desctype == 3) then
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic) * exp(-interatomic_distances(i_atomic,j_atomic))
              else if(this%desctype < 0) then
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = -real(this%desctype) * interatomic_vectors(i_atomic,j_atomic,:) * interatomic_distances(i_atomic,j_atomic)**(this%desctype-2)
              else
                 RAISE_ERROR("molecule_lo_d_calc: not implemented descriptortype",error)
              end if
              descriptor_out%x(i_desc)%grad_data(k,:,j_atomic) = -descriptor_out%x(i_desc)%grad_data(k,:,i_atomic)        ! descriptor wrt j_atomic 
!   call print(descriptor_out%x(i_desc)%grad_data(k,:,:))
            end do
       
         endif

      enddo
      deallocate(shifts)
      deallocate(dist_vec)
      deallocate(atomic_index)
      deallocate(interatomic_vectors)
      deallocate(interatomic_distances)

      call system_timer('molecule_lo_d_calc')

   endsubroutine molecule_lo_d_calc

   subroutine alex_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(alex), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor

      integer :: i, j, n, d, p, q, r, a, b, c, n_i, n_radial, pp, i_desc, l_n_neighbours, desc_index, n_cross, n_descriptors
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij
      real(dp), dimension(3) :: d_ij
      real(dp), dimension(:), allocatable :: neighbour_dists
      real(dp), dimension(:,:), allocatable :: neighbour_vecs
      integer, dimension(116) :: species_map
      real(dp), allocatable :: S0(:), S1(:,:), S2(:,:,:), S0der(:,:,:), S1der(:,:,:,:), S2der(:,:,:,:,:)

      INIT_ERROR(error)

      call system_timer('alex_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("alex_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='alex_calc args_str')) then
            RAISE_ERROR("alex_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("alex_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = alex_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      n_radial = this%power_max-this%power_min+1

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(1))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo


      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         ! number of neighbours for the current atom within the descriptor cutoff
         l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)
         allocate(neighbour_vecs(3,l_n_neighbours), neighbour_dists(l_n_neighbours))
         allocate(S0(n_radial), S1(3,n_radial), S2(3,3,n_radial))
         if(my_do_grad_descriptor) then
            allocate( &
               S0der(n_radial,l_n_neighbours,3), &
               S1der(3,n_radial,l_n_neighbours,3), &
               S2der(3,3,n_radial,l_n_neighbours,3) )
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, diff=d_ij)
            if( r_ij >= this%cutoff ) cycle
            n_i = n_i + 1
            neighbour_vecs(:,n_i) = d_ij
            neighbour_dists(n_i) = r_ij
         end do

         do p = 1,n_radial
            pp = -(p+this%power_min-1)
            S0(p) = sum(neighbour_dists**pp)
            if(my_do_grad_descriptor) then
               do n_i = 1, l_n_neighbours
                  S0der(p,n_i,:) = pp * neighbour_dists(n_i)**(pp-2) * neighbour_vecs(:,n_i)
               enddo
            endif

            S1(:,p) = matmul(neighbour_vecs, neighbour_dists**pp)
            !do a = 1,3
               !S1(a, p) = sum(neighbour_vecs(a,:)*neighbour_dists**pp)
            !end do
            if(my_do_grad_descriptor) then
               do n_i = 1, l_n_neighbours
                  do a = 1,3
                     S1der(a,p,n_i,:) = pp * neighbour_dists(n_i)**(pp-2) * neighbour_vecs(a,n_i) * neighbour_vecs(:,n_i)
                     S1der(a,p,n_i,a) = S1der(a,p,n_i,a) + neighbour_dists(n_i)**pp
                  end do
               enddo
            endif

            !do a=1,3
            do b=1,3
               S2(:,b,p) = matmul(neighbour_vecs, neighbour_vecs(b,:)*neighbour_dists**pp)
               !S2(a,b,p) = sum(neighbour_vecs(a,:)*neighbour_vecs(b,:)*neighbour_dists**pp)
            end do
            !end do

            if(my_do_grad_descriptor) then
               do n_i = 1, l_n_neighbours
                  do a = 1,3
                     do b = 1,3
                        S2der(a,b,p,n_i,:) = pp * neighbour_dists(n_i)**(pp-2) * neighbour_vecs(a,n_i) * neighbour_vecs(b,n_i) * neighbour_vecs(:,n_i)
                     end do
                  end do

                  do a = 1,3
                     do b = 1,3
                        S2der(a,b,p,n_i,b) = S2der(a,b,p,n_i,b) + neighbour_dists(n_i)**pp * neighbour_vecs(a,n_i)
                        S2der(a,b,p,n_i,a) = S2der(a,b,p,n_i,a) + neighbour_dists(n_i)**pp * neighbour_vecs(b,n_i)
                     end do
                  end do
               enddo
            endif
         end do

         descriptor_out%x(i_desc)%data(1:n_radial) = S0
         descriptor_out%x(i_desc)%data(n_radial+1:n_radial+n_radial**2) = reshape(matmul(transpose(S1), S1), (/n_radial**2/))
         desc_index = n_radial+n_radial**2+1
         do p = 1,n_radial
            do q = 1,n_radial
               descriptor_out%x(i_desc)%data(desc_index) = sum(S2(:,:,p) * S2(:,:,q))
               desc_index = desc_index + 1
            end do
         end do
         
         do p = 1,n_radial
            do q = 1,n_radial
               do r = 1,n_radial
                  descriptor_out%x(i_desc)%data(desc_index) =  dot_product(S1(:,p), matmul(S2(:,:,q), S1(:,r)))
                  desc_index = desc_index + 1
               end do
            end do
         end do

         if(my_do_grad_descriptor) then
            n_i = 0
            do n = 1, n_neighbours(at,i)
               j = neighbour(at, i, n, distance = r_ij, shift=shift_ij)
               if( r_ij >= this%cutoff ) cycle

               n_i = n_i + 1

               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.

               descriptor_out%x(i_desc)%grad_data(1:n_radial,:,n_i) = S0der(:,n_i,:)

               desc_index = n_radial + 1
               do p = 1,n_radial
                  do q = 1,n_radial
                     do a = 1, 3
                        do c = 1, 3
                           descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) = descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) + &
                              S1der(a,p,n_i,c)*S1(a,q) + S1(a,p)*S1der(a,q,n_i,c)
                        enddo
                     enddo
                     desc_index = desc_index + 1
                  enddo
               enddo

               do p = 1, n_radial
                  do q = 1, n_radial
                     do a = 1, 3
                        do b = 1, 3
                           do c = 1, 3
                              descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) = descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) + &
                                 S2der(a,b,p,n_i,c)*S2(a,b,q) + S2(a,b,p)*S2der(a,b,q,n_i,c)
                           enddo
                        enddo
                     enddo
                     desc_index = desc_index + 1
                  enddo
               enddo

               do p = 1, n_radial
                  do q = 1, n_radial
                     do r = 1, n_radial
                        do a = 1, 3
                           do b = 1, 3
                              do c = 1, 3
                                 descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) = descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) + &
                                    S1der(a,p,n_i,c) * S2(a,b,q)        * S1(b,r) + &
                                    S1(a,p)        * S2der(a,b,q,n_i,c) * S1(b,r) + &
                                    S1(a,p)        * S2(a,b,q)        * S1der(b,r,n_i,c)
                              enddo
                           enddo
                        enddo
                        desc_index = desc_index + 1
                     enddo
                  enddo
               enddo
            enddo

            descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n_i)
            deallocate(S0der, S1der, S2der)
         endif

         deallocate(neighbour_vecs, neighbour_dists, S0, S1, S2)
      enddo


      call system_timer('alex_calc')

   endsubroutine alex_calc

   subroutine distance_Nb_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(distance_Nb), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, i_desc, i_data, i, j, ii, jj, kk, ll, iConnectivity
      integer, dimension(3) :: s_i, s_j
      real(dp) :: r_ij, fcut_connectivity
      real(dp), dimension(3) :: dfcut_connectivity
      real(dp), dimension(3) :: d_ij
      integer, dimension(:,:,:), allocatable :: atoms_in_descriptors
      real(dp), dimension(:,:), allocatable :: fcut_pair, dfcut_pair
      real(dp), dimension(:,:,:), allocatable :: directions

      logical, dimension(:), pointer :: atom_mask_pointer => null()

      INIT_ERROR(error)

      call system_timer('distance_Nb_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("distance_Nb_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='distance_Nb_calc args_str')) then
            RAISE_ERROR("distance_Nb_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         atom_mask_pointer => null()

         if( has_atom_mask_name ) then
            if( .not. this%compact_clusters ) then
               RAISE_ERROR("distance_Nb_calc: MPI ready only for compact_clusters=T type of distance_Nb.", error)
            endif

            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("distance_2b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      d = distance_Nb_dimensions(this,error)
      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask = atom_mask_pointer, error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross, error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(this%order))
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,this%order))
            allocate(descriptor_out%x(i)%ii(this%order))
            allocate(descriptor_out%x(i)%pos(3,this%order))
            allocate(descriptor_out%x(i)%has_grad_data(this%order))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,this%order))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      if(associated(atom_mask_pointer)) then
         call distance_Nb_calc_get_clusters(this,at,atoms_in_descriptors=atoms_in_descriptors,mask=atom_mask_pointer,error=error)
      else
         call distance_Nb_calc_get_clusters(this,at,atoms_in_descriptors=atoms_in_descriptors,error=error)
      endif

      allocate(fcut_pair(this%order,this%order))
      if( my_do_grad_descriptor ) then
         allocate(dfcut_pair(this%order,this%order), directions(3,this%order,this%order))
      endif

      do i_desc = 1, n_descriptors
         if( this%order == 1 ) then
            descriptor_out%x(i_desc)%data = 0.0_dp
            if( my_do_grad_descriptor ) descriptor_out%x(i_desc)%grad_data = 0.0_dp
         else
            i_data = 0
            do ii = 1, this%order
               i = atoms_in_descriptors(1,ii,i_desc)
               s_i = atoms_in_descriptors(2:4,ii,i_desc)
               do jj = ii+1, this%order
                  i_data = i_data + 1
                  j = atoms_in_descriptors(1,jj,i_desc)
                  s_j = atoms_in_descriptors(2:4,jj,i_desc)
                  d_ij = at%pos(:,j) - at%pos(:,i) + matmul(at%lattice,s_j-s_i)
                  r_ij = sqrt(sum(d_ij**2))

                  fcut_pair(jj,ii) = coordination_function(r_ij,this%cutoff,this%cutoff_transition_width)
                  fcut_pair(ii,jj) = fcut_pair(jj,ii)

                  descriptor_out%x(i_desc)%data(i_data) = r_ij
                  if( my_do_grad_descriptor ) then
                     dfcut_pair(ii,jj) = dcoordination_function(r_ij,this%cutoff,this%cutoff_transition_width)
                     dfcut_pair(jj,ii) = dfcut_pair(ii,jj)

                     directions(:,ii,jj) = d_ij / r_ij
                     directions(:,jj,ii) =  - directions(:,ii,jj)
                     descriptor_out%x(i_desc)%grad_data(i_data,:,jj) = directions(:,ii,jj)
                     descriptor_out%x(i_desc)%grad_data(i_data,:,ii) = &
                       - descriptor_out%x(i_desc)%grad_data(i_data,:,jj)
                  endif
               enddo
            enddo

            descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp
            if ( this%compact_clusters ) then

               descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp

               do jj = 2, this%order
                  descriptor_out%x(i_desc)%covariance_cutoff = descriptor_out%x(i_desc)%covariance_cutoff * fcut_pair(jj,1)
               enddo

               if( my_do_grad_descriptor ) then
                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = 0.0_dp
                  do kk = 2, this%order
                     descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) = 1.0_dp
                     do jj = 2, this%order
                        if( jj == kk ) then
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) = &
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) * dfcut_pair(jj,1) * (-directions(:,jj,1))
                        else
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) = &
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) * fcut_pair(jj,1)
                        endif
                     enddo
                     descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = &
                     descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) - descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk)
                  enddo
               endif


            else
               do iConnectivity = 1, size(this%monomerConnectivities,3)

                  fcut_connectivity = 1.0_dp

                  do ii = 1, this%order
                     do jj = ii+1, this%order
                        if( this%monomerConnectivities(jj,ii,iConnectivity) ) then
                           fcut_connectivity = fcut_connectivity * fcut_pair(jj,ii)
                        else
                           fcut_connectivity = fcut_connectivity * ( 1.0_dp - fcut_pair(jj,ii) )
                        endif
                     enddo
                  enddo
                  descriptor_out%x(i_desc)%covariance_cutoff = descriptor_out%x(i_desc)%covariance_cutoff + fcut_connectivity

                  if( my_do_grad_descriptor ) then
                     do kk = 1, this%order
                        do ll = kk+1, this%order
                           dfcut_connectivity = 1.0_dp
                           do ii = 1, this%order
                              do jj = ii+1, this%order
                                 if( this%monomerConnectivities(jj,ii,iConnectivity) ) then
                                    if( kk == ii .and. ll == jj ) then
                                       dfcut_connectivity = dfcut_connectivity * dfcut_pair(jj,ii) * directions(:,ll,kk)
                                    elseif( kk == jj .and. ll == ii ) then
                                       dfcut_connectivity = dfcut_connectivity * dfcut_pair(jj,ii) * directions(:,ll,kk)
                                    else
                                       dfcut_connectivity = dfcut_connectivity * fcut_pair(jj,ii)
                                    endif
                                 else
                                    if( kk == ii .and. ll == jj ) then
                                       dfcut_connectivity = - dfcut_connectivity * dfcut_pair(jj,ii) * directions(:,ll,kk)
                                    elseif( kk == jj .and. ll == ii) then
                                       dfcut_connectivity = - dfcut_connectivity * dfcut_pair(jj,ii) * directions(:,ll,kk)
                                    else
                                       dfcut_connectivity = dfcut_connectivity * ( 1.0_dp - fcut_pair(jj,ii) )
                                    endif
                                 endif
                              enddo !jj
                           enddo !ii
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) + &
                              dfcut_connectivity
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,ll) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,ll) - &
                              dfcut_connectivity
                        enddo !ll
                     enddo !kk
                  endif

               enddo
            endif

         endif

         descriptor_out%x(i_desc)%ci = atoms_in_descriptors(1,:,i_desc)
         descriptor_out%x(i_desc)%has_data = .true.
         if( my_do_grad_descriptor ) then
            descriptor_out%x(i_desc)%ii = descriptor_out%x(i_desc)%ci
            descriptor_out%x(i_desc)%pos = at%pos(:,descriptor_out%x(i_desc)%ii) + &
               matmul(at%lattice,atoms_in_descriptors(2:4,:,i_desc))
            descriptor_out%x(i_desc)%has_grad_data = .true.
         endif

      enddo

      if(allocated(atoms_in_descriptors)) deallocate(atoms_in_descriptors)
      if(allocated(fcut_pair)) deallocate(fcut_pair)
      if(allocated(dfcut_pair)) deallocate(dfcut_pair)
      if(allocated(directions)) deallocate(directions)

      call system_timer('distance_Nb_calc')

   endsubroutine distance_Nb_calc

   subroutine distance_Nb_calc_get_clusters(this,at,atoms_in_descriptors,n_descriptors,mask,error)

      type(distance_Nb), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, dimension(:,:,:), intent(out), allocatable, optional :: atoms_in_descriptors
      integer, intent(out), optional :: n_descriptors
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: error

      integer, dimension(:,:,:), allocatable :: my_atoms_in_descriptors

      if( present(atoms_in_descriptors) ) then
         call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors,n_descriptors=n_descriptors,mask=mask,error=error)
      else
         call distance_Nb_calc_neighbour_loop(this,at,my_atoms_in_descriptors,n_descriptors=n_descriptors,mask=mask,error=error)
         if(allocated(my_atoms_in_descriptors)) deallocate(my_atoms_in_descriptors)
      endif

   endsubroutine distance_Nb_calc_get_clusters

!   recursive subroutine distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors,n_descriptors,error)
!
!      type(distance_Nb), intent(in) :: this
!      type(atoms), intent(in) :: at
!      integer, dimension(:,:,:), intent(inout), allocatable :: atoms_in_descriptors
!      integer, intent(out), optional :: n_descriptors
!      integer, intent(out), optional :: error
!
!      integer, save :: current_order = 0
!      integer :: i, j, n, order, i_desc, d
!      real(dp) :: r_ij
!      integer, dimension(3) :: shift_i, shift_j, shift_ij
!      integer, dimension(:,:), allocatable :: current_descriptor
!
!      type(LinkedList_i2d), pointer :: LL_atoms_in_descriptors => null()
!
!      INIT_ERROR(error)
!
!      current_order = current_order + 1
!
!      if( current_order == 1 ) then
!         allocate(current_descriptor(4,1))
!
!         do i = 1, at%N
!            if( any( at%Z(i) == this%Z ) .or. any( 0 == this%Z ) ) then
!               current_descriptor(:,1) = (/i,0,0,0/)
!               call append(LL_atoms_in_descriptors,current_descriptor,error)
!            endif
!         enddo
!
!         deallocate(current_descriptor)
!         call retrieve(LL_atoms_in_descriptors,atoms_in_descriptors)
!         call finalise(LL_atoms_in_descriptors)
!         if( this%order > 1 ) &
!            call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors = atoms_in_descriptors,n_descriptors=n_descriptors,error=error)
!
!         if( present(n_descriptors) ) n_descriptors = size(atoms_in_descriptors,3)
!      else
!         if( .not. allocated(atoms_in_descriptors) ) then
!            RAISE_ERROR("distance_Nb_calc_neighbour_loop: atoms_in_descriptors must be allocated",error)
!         endif
!
!         allocate(current_descriptor(4,current_order))
!         do i_desc = 1, size(atoms_in_descriptors,3)
!            do order = 1, size(atoms_in_descriptors,2)
!               i = atoms_in_descriptors(1,order,i_desc)
!               shift_i = atoms_in_descriptors(2:4,order,i_desc)
!               loop_n: do n = 1, n_neighbours(at,i)
!                  j = neighbour(at,i,n,distance = r_ij, shift = shift_ij)
!
!                  if( r_ij > this%cutoff ) cycle
!                  if( .not. is_subset(this%Z, at%Z( (/j,atoms_in_descriptors(1,:,i_desc)/) ), error) .and. all(this%Z /= 0) ) cycle
!
!                  shift_j = shift_ij + shift_i
!
!                  current_descriptor(:,1:current_order-1) = atoms_in_descriptors(:,:,i_desc)
!                  current_descriptor(:,current_order) = (/j, shift_j/)
!                  if( order_and_check_for_duplicates(current_descriptor,at) ) then
!                     do d = current_order, 1, -1
!                        current_descriptor(2:4,d) = current_descriptor(2:4,d) - current_descriptor(2:4,1)
!                     enddo
!                     if( .not. is_in_LinkedList(LL_atoms_in_descriptors,current_descriptor,error) ) &
!                       call append(LL_atoms_in_descriptors,current_descriptor,error)
!                  endif
!               enddo loop_n
!            enddo
!         enddo
!
!         deallocate(current_descriptor)
!         call retrieve(LL_atoms_in_descriptors,atoms_in_descriptors)
!         call finalise(LL_atoms_in_descriptors)
!         if( current_order < this%order ) &
!            call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors = atoms_in_descriptors,n_descriptors=n_descriptors,error=error)
!      endif
!
!      current_order = current_order - 1
!
!   endsubroutine distance_Nb_calc_neighbour_loop

   recursive subroutine distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors,n_descriptors,mask,error)

      type(distance_Nb), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, dimension(:,:,:), intent(inout), allocatable :: atoms_in_descriptors
      integer, intent(out), optional :: n_descriptors
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: error

      integer, save :: current_order = 0
      integer :: i, j, n, order, i_desc, d
      real(dp) :: r_ij
      integer, dimension(3) :: shift_i, shift_j, shift_ij
      integer, dimension(:,:), allocatable :: current_descriptor

      type(Table)  :: Table_atoms_in_descriptors, Table_atoms_in_descriptors_uniq

      INIT_ERROR(error)

      current_order = current_order + 1

      if( current_order == 1 ) then
         call initialise(Table_atoms_in_descriptors, Nint = 4*current_order, Nreal = 0, Nstr = 0, Nlogical = 0, error=error)
         allocate(current_descriptor(4,1))

         do i = 1, at%N
            if( any( at%Z(i) == this%Z ) .or. any( 0 == this%Z ) ) then
               if( present(mask) ) then
                  if( .not. mask(i) ) cycle
               endif

               current_descriptor(:,1) = (/i,0,0,0/)
               call append(Table_atoms_in_descriptors,current_descriptor(:,1))
            endif
         enddo

         deallocate(current_descriptor)

         allocate(atoms_in_descriptors(4,1,Table_atoms_in_descriptors%N))
         atoms_in_descriptors = reshape(Table_atoms_in_descriptors%int(:,1:Table_atoms_in_descriptors%N),(/4,1,Table_atoms_in_descriptors%N/))
         
         call finalise(Table_atoms_in_descriptors)

         if( this%order > 1 ) &
            call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors = atoms_in_descriptors,n_descriptors=n_descriptors,error=error)

         if( present(n_descriptors) ) n_descriptors = size(atoms_in_descriptors,3)
      else
         if( .not. allocated(atoms_in_descriptors) ) then
            RAISE_ERROR("distance_Nb_calc_neighbour_loop: atoms_in_descriptors must be allocated",error)
         endif

         call initialise(Table_atoms_in_descriptors, Nint = 4*current_order, Nreal = 0, Nstr = 0, Nlogical = 0, error=error)
         allocate(current_descriptor(4,current_order))

         do i_desc = 1, size(atoms_in_descriptors,3)
            do order = 1, merge(1,size(atoms_in_descriptors,2),this%compact_clusters) !size(atoms_in_descriptors,2)
               ! if compact_clusters == T, only neighbours of the first (central) atom is considered
               i = atoms_in_descriptors(1,order,i_desc)
               shift_i = atoms_in_descriptors(2:4,order,i_desc)
               loop_n: do n = 1, n_neighbours(at,i)
                  j = neighbour(at,i,n,distance = r_ij, shift = shift_ij)

                  if( r_ij > this%cutoff ) cycle
                  if( .not. is_subset(this%Z, at%Z( (/j,atoms_in_descriptors(1,:,i_desc)/) ), error) .and. all(this%Z /= 0) ) cycle

                  shift_j = shift_ij + shift_i

                  current_descriptor(:,1:current_order-1) = atoms_in_descriptors(:,:,i_desc)
                  current_descriptor(:,current_order) = (/j, shift_j/)
                  if( order_and_check_for_duplicates(current_descriptor(:,merge(2,1,this%compact_clusters):),at) ) then
                     ! if compact_clusters == T, leave first atom alone
                     do d = current_order, 1, -1
                        current_descriptor(2:4,d) = current_descriptor(2:4,d) - current_descriptor(2:4,1)
                     enddo
                     call append(Table_atoms_in_descriptors,reshape(current_descriptor,(/4*current_order/)))

                     !if( .not. is_in_LinkedList(LL_atoms_in_descriptors,current_descriptor,error) ) &
                     !  call append(LL_atoms_in_descriptors,current_descriptor,error)
                  endif
               enddo loop_n
            enddo
         enddo

         deallocate(current_descriptor,atoms_in_descriptors)
         call initialise(Table_atoms_in_descriptors_uniq, Nint = 4*current_order, Nreal = 0, Nstr = 0, Nlogical = 0, error=error)

         if( Table_atoms_in_descriptors%N > 0 ) then
            call heap_sort(Table_atoms_in_descriptors%int(:,1:Table_atoms_in_descriptors%N))
            call append(Table_atoms_in_descriptors_uniq,Table_atoms_in_descriptors%int(:,1))
            do i_desc = 2, Table_atoms_in_descriptors%N
               if( .not. all( Table_atoms_in_descriptors%int(:,i_desc) == Table_atoms_in_descriptors%int(:,i_desc-1) ) ) &
                   call append(Table_atoms_in_descriptors_uniq,Table_atoms_in_descriptors%int(:,i_desc))
            enddo
         endif

         allocate(atoms_in_descriptors(4,current_order,Table_atoms_in_descriptors_uniq%N))
         atoms_in_descriptors = reshape(Table_atoms_in_descriptors_uniq%int(:,1:Table_atoms_in_descriptors_uniq%N),(/4,current_order,Table_atoms_in_descriptors_uniq%N/))

         call finalise(Table_atoms_in_descriptors)
         call finalise(Table_atoms_in_descriptors_uniq)

         if( current_order < this%order ) &
            call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors = atoms_in_descriptors,n_descriptors=n_descriptors,error=error)
      endif

      current_order = current_order - 1

   endsubroutine distance_Nb_calc_neighbour_loop

   function order_and_check_for_duplicates(array,at)
      integer, dimension(:,:), intent(inout) :: array
      type(atoms), intent(in) :: at
      logical :: order_and_check_for_duplicates

      integer :: ii, jj, n
      integer, dimension(size(array,1)) :: tmp
      logical :: do_swap

      integer, dimension(size(array,1)+1,size(array,2)) :: Z_array

      Z_array(1,:) = at%Z(array(1,:))
      Z_array(2:,:) = array(:,:)

      call heap_sort(Z_array)

      do ii = 2, size(Z_array,2)
         if( all( Z_array(:,ii-1) == Z_array(:,ii) ) ) then
            order_and_check_for_duplicates = .false.
            return
         endif
      enddo

      array(:,:) = Z_array(2:,:)

      order_and_check_for_duplicates = .true.

   endfunction order_and_check_for_duplicates

   function is_subset(set,subset,error)
      logical :: is_subset
      integer, dimension(:), intent(in) :: set, subset
      integer, optional, intent(out) :: error

      logical, dimension(size(set)) :: found
      integer :: i, j

      INIT_ERROR(error)
      if( size(set) < size(subset) ) then
         RAISE_ERROR("is_subset: size of set must be greater than or equal to the size of subset",error)
      endif

      found = .false.
      loop_i: do i = 1, size(subset)
         do j = 1, size(set)
            if(set(j) == subset(i) .and. .not. found(j)) then
               found(j) = .true.
               cycle loop_i
            endif
         enddo
      enddo loop_i

      is_subset = ( count(found) == size(subset) )

   endfunction is_subset

   function descriptor_dimensions(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: descriptor_dimensions

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            descriptor_dimensions = bispectrum_SO4_dimensions(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            descriptor_dimensions = bispectrum_SO3_dimensions(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            descriptor_dimensions = behler_dimensions(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            descriptor_dimensions = distance_2b_dimensions(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            descriptor_dimensions = coordination_dimensions(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            descriptor_dimensions = angle_3b_dimensions(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            descriptor_dimensions = co_angle_3b_dimensions(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            descriptor_dimensions = co_distance_2b_dimensions(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            descriptor_dimensions = cosnx_dimensions(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            descriptor_dimensions = trihis_dimensions(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            descriptor_dimensions = water_monomer_dimensions(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            descriptor_dimensions = water_dimer_dimensions(this%descriptor_water_dimer,error)
         case(DT_A2_DIMER)
            descriptor_dimensions = A2_dimer_dimensions(this%descriptor_A2_dimer,error)
         case(DT_AB_DIMER)
            descriptor_dimensions = AB_dimer_dimensions(this%descriptor_AB_dimer,error)
         case(DT_BOND_REAL_SPACE)
            descriptor_dimensions = bond_real_space_dimensions(this%descriptor_bond_real_space,error)
         case(DT_ATOM_REAL_SPACE)
            descriptor_dimensions = atom_real_space_dimensions(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            descriptor_dimensions = power_so3_dimensions(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            descriptor_dimensions = power_so4_dimensions(this%descriptor_power_so4,error)
         case(DT_SOAP)
            descriptor_dimensions = soap_dimensions(this%descriptor_soap,error)
         case(DT_AN_MONOMER)
            descriptor_dimensions = AN_monomer_dimensions(this%descriptor_AN_monomer,error)
         case(DT_GENERAL_MONOMER)
            descriptor_dimensions = general_monomer_dimensions(this%descriptor_general_monomer,error)
         case(DT_GENERAL_DIMER)
            descriptor_dimensions = general_dimer_dimensions(this%descriptor_general_dimer,error)
         case(DT_GENERAL_TRIMER)
            descriptor_dimensions = general_trimer_dimensions(this%descriptor_general_trimer,error)
         case(DT_RDF)
            descriptor_dimensions = rdf_dimensions(this%descriptor_rdf,error)
         case(DT_AS_DISTANCE_2b)
            descriptor_dimensions = as_distance_2b_dimensions(this%descriptor_as_distance_2b,error)
         case(DT_MOLECULE_LO_D)
            descriptor_dimensions = molecule_lo_d_dimensions(this%descriptor_molecule_lo_d,error)
         case(DT_ALEX)
            descriptor_dimensions = alex_dimensions(this%descriptor_alex,error)
         case(DT_COM_DIMER)
            descriptor_dimensions = com_dimer_dimensions(this%descriptor_com_dimer,error)
         case(DT_DISTANCE_Nb)
            descriptor_dimensions = distance_Nb_dimensions(this%descriptor_distance_Nb,error)
         case default
            RAISE_ERROR("descriptor_dimensions: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_dimensions

   function bispectrum_SO4_dimensions(this,error) result(i)
      type(bispectrum_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i
      integer :: j, j1, j2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_dimensions: descriptor object not initialised", error)
      endif

      i = 0
      do j1 = 0, this%j_max
         j2 = j1
         !do j2 = 0, this%j_max
            do j = abs(j1-j2), min(this%j_max,j1+j2)
               if( mod(j1+j2+j,2) == 1 ) cycle
               i = i + 1
            enddo
         !enddo
      enddo

   endfunction bispectrum_SO4_dimensions

   function bispectrum_SO3_dimensions(this,error) result(i)
      type(bispectrum_SO3), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i
      integer :: a, l1, l2, l

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_dimensions: descriptor object not initialised", error)
      endif

      i = 0
      do a = 1, this%n_max
         do l1 = 0, this%l_max
            l2 = l1
            !do l2 = 0, this%l_max
               do l = abs(l1-l2), min(this%l_max,l1+l2)
                  if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                  i = i + 1
               enddo
            !enddo
         enddo
      enddo

   endfunction bispectrum_SO3_dimensions

   function behler_dimensions(this,error) result(i)
      type(behler), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_g2 + this%n_g3

   endfunction behler_dimensions

   function distance_2b_dimensions(this,error) result(i)
      type(distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_dimensions: descriptor object not initialised", error)
      endif

      i = 1

   endfunction distance_2b_dimensions

   function coordination_dimensions(this,error) result(i)
      type(coordination), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_dimensions: descriptor object not initialised", error)
      endif

      i = 1

   endfunction coordination_dimensions

   function angle_3b_dimensions(this,error) result(i)
      type(angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction angle_3b_dimensions

   function co_angle_3b_dimensions(this,error) result(i)
      type(co_angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_dimensions: descriptor object not initialised", error)
      endif

      i = 4

   endfunction co_angle_3b_dimensions

   function co_distance_2b_dimensions(this,error) result(i)
      type(co_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction co_distance_2b_dimensions

   function cosnx_dimensions(this,error) result(i)
      type(cosnx), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_max*(this%l_max+1)

   endfunction cosnx_dimensions

   function trihis_dimensions(this,error) result(i)
      type(trihis), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_gauss

   endfunction trihis_dimensions

   function water_monomer_dimensions(this,error) result(i)
      type(water_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction water_monomer_dimensions

   function water_dimer_dimensions(this,error) result(i)
      type(water_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 15

   endfunction water_dimer_dimensions

   function A2_dimer_dimensions(this,error) result(i)
      type(A2_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 6

   endfunction A2_dimer_dimensions

   function AB_dimer_dimensions(this,error) result(i)
      type(AB_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 6

   endfunction AB_dimer_dimensions

   function bond_real_space_dimensions(this,error) result(i)
      type(bond_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_dimensions: descriptor object not initialised", error)
      endif

      i = 2 + (1 + 2 * this%max_neighbours) * this%max_neighbours

   endfunction bond_real_space_dimensions

   function atom_real_space_dimensions(this,error) result(i)
      type(atom_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_dimensions: descriptor object not initialised", error)
      endif

      i = 2 * (this%l_max+1)**2 + 2

   endfunction atom_real_space_dimensions

   function power_so3_dimensions(this,error) result(i)
      type(power_so3), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_max*(this%l_max+1)

   endfunction power_so3_dimensions

   function power_SO4_dimensions(this,error) result(i)
      type(power_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_dimensions: descriptor object not initialised", error)
      endif

      i = this%j_max + 1

   endfunction power_SO4_dimensions

   function soap_dimensions(this,error) result(i)
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_dimensions: descriptor object not initialised", error)
      endif

      !i = (this%l_max+1) * ( this%n_max * (this%n_max+1) / 2 ) * ( this%n_species * (this%n_species+1) / 2 ) + 1
      !i = (this%l_max+1) * this%n_max**2 * this%n_species**2 + 1
      if(this%diagonal_radial) then
         i = (this%l_max+1) * this%n_max * this%n_species * (this%n_species+1) / 2 + 1
      else
         i = (this%l_max+1) * ( (this%n_max*this%n_species)*(this%n_max*this%n_species+1) ) / 2 + 1
      endif

   endfunction soap_dimensions

   function AN_monomer_dimensions(this,error) result(i)
      type(AN_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AN_monomer_dimensions: descriptor object not initialised", error)
      endif

      i = this%N * (this%N - 1) / 2

   endfunction AN_monomer_dimensions

   function general_monomer_dimensions(this,error) result(i)
      type(general_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_monomer_dimensions: descriptor object not initialised", error)
      endif
      if(.not. this%permutation_data%initialised) then
         RAISE_ERROR("general_monomer_dimensions: descriptor object's permutation data not initialised", error)
      endif

      i = size(this%permutation_data%dist_vec)

   endfunction general_monomer_dimensions

   function com_dimer_dimensions(this,error) result(i)
      type(com_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("com_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 1

   endfunction com_dimer_dimensions

   function general_dimer_dimensions(this,error) result(i)
      type(general_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_dimer_dimensions: descriptor object not initialised", error)
      endif
      if(.not. this%permutation_data%initialised) then
         RAISE_ERROR("general_dimer_dimensions: descriptor object's permutation data not initialised", error)
      endif

      i = size(this%permutation_data%dist_vec)

   endfunction general_dimer_dimensions


   function general_trimer_dimensions(this,error) result(i)
      type(general_trimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_trimer_dimensions: descriptor object not initialised", error)
      endif
      if(.not. this%permutation_data%initialised) then
         RAISE_ERROR("general_trimer_dimensions: descriptor object's permutation data not initialised", error)
      endif

      i = size(this%permutation_data%dist_vec)

   endfunction general_trimer_dimensions

   function rdf_dimensions(this,error) result(i)
      type(rdf), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("rdf_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_gauss

   endfunction rdf_dimensions

   function as_distance_2b_dimensions(this,error) result(i)
      type(as_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("as_distance_2b_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction as_distance_2b_dimensions

   function molecule_lo_d_dimensions(this,error) result(i)
      type(molecule_lo_d), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("molecule_lo_d_dimensions: descriptor object not initialised", error)
      endif
      if(.not. this%permutation_data%initialised) then
         RAISE_ERROR("molecule_lo_d_dimensions: descriptor object's permutation data not initialised", error)
      endif

      i = size(this%included_components)

   endfunction molecule_lo_d_dimensions

   function alex_dimensions(this,error) result(i)
      type(alex), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i, nradial

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("alex_dimensions: descriptor object not initialised", error)
      endif

      nradial = this%power_max-this%power_min + 1
      i = nradial+2*nradial**2+nradial**3

   endfunction alex_dimensions

   function distance_Nb_dimensions(this,error) result(i)
      type(distance_Nb), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_Nb_dimensions: descriptor object not initialised", error)
      endif

      i = max(1,this%order * ( this%order - 1 ) / 2)

   endfunction distance_Nb_dimensions

   function descriptor_cutoff(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: descriptor_cutoff

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            descriptor_cutoff = cutoff(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            descriptor_cutoff = cutoff(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            descriptor_cutoff = cutoff(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            descriptor_cutoff = cutoff(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            descriptor_cutoff = cutoff(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            descriptor_cutoff = cutoff(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            descriptor_cutoff = cutoff(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            descriptor_cutoff = cutoff(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            descriptor_cutoff = cutoff(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            descriptor_cutoff = cutoff(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            descriptor_cutoff = cutoff(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_water_dimer,error)
         case(DT_A2_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_A2_dimer,error)
         case(DT_AB_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_AB_dimer,error)
         case(DT_BOND_REAL_SPACE)
            descriptor_cutoff = cutoff(this%descriptor_bond_real_space,error)
         case(DT_ATOM_REAL_SPACE)
            descriptor_cutoff = cutoff(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            descriptor_cutoff = cutoff(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            descriptor_cutoff = cutoff(this%descriptor_power_so4,error)
         case(DT_SOAP)
            descriptor_cutoff = cutoff(this%descriptor_soap,error)
         case(DT_AN_MONOMER)
            descriptor_cutoff = cutoff(this%descriptor_AN_monomer,error)
         case(DT_GENERAL_MONOMER)
            descriptor_cutoff = cutoff(this%descriptor_general_monomer,error)
         case(DT_GENERAL_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_general_dimer,error)
         case(DT_GENERAL_TRIMER)
            descriptor_cutoff = cutoff(this%descriptor_general_trimer,error)
         case(DT_RDF)
            descriptor_cutoff = cutoff(this%descriptor_rdf,error)
         case(DT_MOLECULE_LO_D)
            descriptor_cutoff = cutoff(this%descriptor_molecule_lo_d,error)
         case(DT_ALEX)
            descriptor_cutoff = cutoff(this%descriptor_alex,error)
         case(DT_COM_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_com_dimer,error)
         case(DT_DISTANCE_Nb)
            descriptor_cutoff = cutoff(this%descriptor_distance_Nb,error)
         case default
            RAISE_ERROR("descriptor_cutoff: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_cutoff

   function bispectrum_SO4_cutoff(this,error) 
      type(bispectrum_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bispectrum_SO4_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_cutoff: descriptor object not initialised", error)
      endif

      bispectrum_SO4_cutoff = this%cutoff

   endfunction bispectrum_SO4_cutoff

   function bispectrum_SO3_cutoff(this,error) 
      type(bispectrum_SO3), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bispectrum_SO3_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_cutoff: descriptor object not initialised", error)
      endif

      bispectrum_SO3_cutoff = this%cutoff

   endfunction bispectrum_SO3_cutoff

   function behler_cutoff(this,error) 
      type(behler), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: behler_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_cutoff: descriptor object not initialised", error)
      endif

      behler_cutoff = this%cutoff

   endfunction behler_cutoff

   function distance_2b_cutoff(this,error) 
      type(distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: distance_2b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_cutoff: descriptor object not initialised", error)
      endif

      distance_2b_cutoff = this%cutoff

   endfunction distance_2b_cutoff

   function co_distance_2b_cutoff(this,error) 
      type(co_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: co_distance_2b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_cutoff: descriptor object not initialised", error)
      endif

      co_distance_2b_cutoff = this%cutoff

   endfunction co_distance_2b_cutoff

   function coordination_cutoff(this,error) 
      type(coordination), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: coordination_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_cutoff: descriptor object not initialised", error)
      endif

      coordination_cutoff = this%cutoff

   endfunction coordination_cutoff

   function angle_3b_cutoff(this,error) 
      type(angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: angle_3b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_cutoff: descriptor object not initialised", error)
      endif

      angle_3b_cutoff = this%cutoff

   endfunction angle_3b_cutoff

   function co_angle_3b_cutoff(this,error) 
      type(co_angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: co_angle_3b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_cutoff: descriptor object not initialised", error)
      endif

      co_angle_3b_cutoff = this%cutoff

   endfunction co_angle_3b_cutoff

   function cosnx_cutoff(this,error) 
      type(cosnx), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: cosnx_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_cutoff: descriptor object not initialised", error)
      endif

      cosnx_cutoff = this%cutoff

   endfunction cosnx_cutoff

   function trihis_cutoff(this,error) 
      type(trihis), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: trihis_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_cutoff: descriptor object not initialised", error)
      endif

      trihis_cutoff = this%cutoff

   endfunction trihis_cutoff

   function water_monomer_cutoff(this,error) 
      type(water_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: water_monomer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_cutoff: descriptor object not initialised", error)
      endif

      water_monomer_cutoff = this%cutoff

   endfunction water_monomer_cutoff

   function water_dimer_cutoff(this,error) 
      type(water_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: water_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_cutoff: descriptor object not initialised", error)
      endif

      water_dimer_cutoff = this%cutoff

   endfunction water_dimer_cutoff

   function A2_dimer_cutoff(this,error) 
      type(A2_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: A2_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_cutoff: descriptor object not initialised", error)
      endif

      A2_dimer_cutoff = this%cutoff

   endfunction A2_dimer_cutoff

   function AB_dimer_cutoff(this,error) 
      type(AB_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: AB_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_cutoff: descriptor object not initialised", error)
      endif

      AB_dimer_cutoff = this%cutoff

   endfunction AB_dimer_cutoff

   function bond_real_space_cutoff(this,error)
      type(bond_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bond_real_space_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_cutoff: descriptor object not initialised", error)
      endif

      bond_real_space_cutoff = max(this%cutoff, this%bond_cutoff)

   endfunction bond_real_space_cutoff

   function atom_real_space_cutoff(this,error)
      type(atom_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: atom_real_space_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_cutoff: descriptor object not initialised", error)
      endif

      atom_real_space_cutoff = this%cutoff

   endfunction atom_real_space_cutoff

   function power_so3_cutoff(this,error) 
      type(power_so3), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: power_so3_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_cutoff: descriptor object not initialised", error)
      endif

      power_so3_cutoff = this%cutoff

   endfunction power_so3_cutoff

   function power_so4_cutoff(this,error) 
      type(power_so4), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: power_so4_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so4_cutoff: descriptor object not initialised", error)
      endif

      power_so4_cutoff = this%cutoff

   endfunction power_so4_cutoff

   function soap_cutoff(this,error) 
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: soap_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_cutoff: descriptor object not initialised", error)
      endif

      soap_cutoff = this%cutoff

   endfunction soap_cutoff

   function AN_monomer_cutoff(this,error) 
      type(AN_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: AN_monomer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AN_monomer_cutoff: descriptor object not initialised", error)
      endif

      AN_monomer_cutoff = this%cutoff

   endfunction AN_monomer_cutoff

   function general_monomer_cutoff(this,error) 
      type(general_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: general_monomer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_monomer_cutoff: descriptor object not initialised", error)
      endif

      general_monomer_cutoff = this%cutoff

   endfunction general_monomer_cutoff

   function com_dimer_cutoff(this,error) 
      type(com_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: com_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("com_dimer_cutoff: descriptor object not initialised", error)
      endif

      com_dimer_cutoff = this%cutoff

   endfunction com_dimer_cutoff

   function general_dimer_cutoff(this,error) 
      type(general_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: general_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_dimer_cutoff: descriptor object not initialised", error)
      endif

      general_dimer_cutoff = this%cutoff

   endfunction general_dimer_cutoff

   function general_trimer_cutoff(this,error) 
      type(general_trimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: general_trimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_trimer_cutoff: descriptor object not initialised", error)
      endif

      general_trimer_cutoff = this%cutoff

   endfunction general_trimer_cutoff

   function rdf_cutoff(this,error) 
      type(rdf), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: rdf_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("rdf_cutoff: descriptor object not initialised", error)
      endif

      rdf_cutoff = this%cutoff

   endfunction rdf_cutoff

   function as_distance_2b_cutoff(this,error) 
      type(as_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: as_distance_2b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("as_distance_2b_cutoff: descriptor object not initialised", error)
      endif

      as_distance_2b_cutoff = this%max_cutoff

   endfunction as_distance_2b_cutoff

   function molecule_lo_d_cutoff(this,error) 
      type(molecule_lo_d), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: molecule_lo_d_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("molecule_lo_d_cutoff: descriptor object not initialised", error)
      endif

      molecule_lo_d_cutoff = this%cutoff

   endfunction molecule_lo_d_cutoff

   function alex_cutoff(this,error) 
      type(alex), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: alex_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("alex_cutoff: descriptor object not initialised", error)
      endif

      alex_cutoff = this%cutoff

   endfunction alex_cutoff

   function distance_Nb_cutoff(this,error) 
      type(distance_Nb), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: distance_Nb_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_Nb_cutoff: descriptor object not initialised", error)
      endif

      distance_Nb_cutoff = this%cutoff

   endfunction distance_Nb_cutoff

   subroutine descriptor_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call bispectrum_SO4_sizes(this%descriptor_bispectrum_SO4,at,n_descriptors,n_cross,mask,error=error)
         case(DT_BISPECTRUM_SO3)
            call bispectrum_SO3_sizes(this%descriptor_bispectrum_SO3,at,n_descriptors,n_cross,mask,error=error)
         case(DT_BEHLER)
            call behler_sizes(this%descriptor_behler,at,n_descriptors,n_cross,mask,error=error)
         case(DT_DISTANCE_2b)
            call distance_2b_sizes(this%descriptor_distance_2b,at,n_descriptors,n_cross,mask,error=error)
         case(DT_COORDINATION)
            call coordination_sizes(this%descriptor_coordination,at,n_descriptors,n_cross,mask,error=error)
         case(DT_ANGLE_3B)
            call angle_3b_sizes(this%descriptor_angle_3b,at,n_descriptors,n_cross,mask,error=error)
         case(DT_CO_ANGLE_3B)
            call co_angle_3b_sizes(this%descriptor_co_angle_3b,at,n_descriptors,n_cross,mask,error=error)
         case(DT_CO_DISTANCE_2b)
            call co_distance_2b_sizes(this%descriptor_co_distance_2b,at,n_descriptors,n_cross,mask,error=error)
         case(DT_COSNX)
            call cosnx_sizes(this%descriptor_cosnx,at,n_descriptors,n_cross,mask,error=error)
         case(DT_TRIHIS)
            call trihis_sizes(this%descriptor_trihis,at,n_descriptors,n_cross,mask,error=error)
         case(DT_WATER_MONOMER)
            call water_monomer_sizes(this%descriptor_water_monomer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_WATER_DIMER)
            call water_dimer_sizes(this%descriptor_water_dimer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_A2_DIMER)
            call A2_dimer_sizes(this%descriptor_A2_dimer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_AB_DIMER)
            call AB_dimer_sizes(this%descriptor_AB_dimer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_BOND_REAL_SPACE)
            call bond_real_space_sizes(this%descriptor_bond_real_space,at,n_descriptors,n_cross,mask,error=error)
         case(DT_ATOM_REAL_SPACE)
            call atom_real_space_sizes(this%descriptor_atom_real_space,at,n_descriptors,n_cross,mask,error=error)
         case(DT_POWER_SO3)
            call power_so3_sizes(this%descriptor_power_so3,at,n_descriptors,n_cross,mask,error=error)
         case(DT_POWER_SO4)
            call power_so4_sizes(this%descriptor_power_so4,at,n_descriptors,n_cross,mask,error=error)
         case(DT_SOAP)
            call soap_sizes(this%descriptor_soap,at,n_descriptors,n_cross,mask,error=error)
         case(DT_AN_MONOMER)
            call AN_monomer_sizes(this%descriptor_AN_monomer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_GENERAL_MONOMER)
            call general_monomer_sizes(this%descriptor_general_monomer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_GENERAL_DIMER)
            call general_dimer_sizes(this%descriptor_general_dimer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_GENERAL_TRIMER)
            call general_trimer_sizes(this%descriptor_general_trimer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_RDF)
            call rdf_sizes(this%descriptor_rdf,at,n_descriptors,n_cross,mask,error=error)
         case(DT_AS_DISTANCE_2b)
            call as_distance_2b_sizes(this%descriptor_as_distance_2b,at,n_descriptors,n_cross,mask,error=error)
         case(DT_MOLECULE_LO_D)
            call molecule_lo_d_sizes(this%descriptor_molecule_lo_d,at,n_descriptors,n_cross,mask,error=error)
         case(DT_ALEX)
            call alex_sizes(this%descriptor_alex,at,n_descriptors,n_cross,mask,error=error)
         case(DT_COM_DIMER)
            call com_dimer_sizes(this%descriptor_com_dimer,at,n_descriptors,n_cross,mask,error=error)
         case(DT_DISTANCE_Nb)
            call distance_Nb_sizes(this%descriptor_distance_Nb,at,n_descriptors,n_cross,mask,error=error)
         case default
            RAISE_ERROR("descriptor_sizes: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_sizes
   subroutine bispectrum_SO4_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(bispectrum_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine bispectrum_SO4_sizes

   subroutine bispectrum_SO3_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(bispectrum_SO3), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine bispectrum_SO3_sizes

   subroutine behler_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(behler), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo
   endsubroutine behler_sizes

   subroutine distance_2b_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2
      real(dp) :: r_ij

      logical :: needs_resid
      integer, dimension(:), pointer :: resid_pointer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_sizes: descriptor object not initialised", error)
      endif

      needs_resid = this%only_intra .or. this%only_inter
      if (needs_resid) then
         if (.not. assign_pointer(at, trim(this%resid_name), resid_pointer)) then
            RAISE_ERROR("distance_2b_sizes did not find "//trim(this%resid_name)//" property (residue id) in the atoms object.", error)
         end if
      else
         resid_pointer => null()
      end if

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif

         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance=r_ij)
            if(r_ij >= this%cutoff) cycle
!if(r_ij < 3.5_dp) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            if (needs_resid) then
               if (this%only_intra .and. resid_pointer(i) /= resid_pointer(j)) cycle
               if (this%only_inter .and. resid_pointer(i) == resid_pointer(j)) cycle
            end if

            n_descriptors = n_descriptors + 1
         enddo
      enddo

      n_cross = n_descriptors*2

   endsubroutine distance_2b_sizes

   subroutine coordination_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(coordination), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine coordination_sizes

   subroutine angle_3b_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i, j, k, n, m
      real(dp) :: r_ij, r_ik
      logical :: Zk1, Zk2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij)
            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, n_neighbours(at,i)
               if( n == m ) cycle

               k = neighbour(at, i, m, distance = r_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)
               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               n_descriptors = n_descriptors + 1
            enddo
         enddo
      enddo
      n_cross = n_descriptors * 3

   endsubroutine angle_3b_sizes

   subroutine co_angle_3b_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(co_angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i, j, k, n, m, n_neighbours_coordination
      real(dp) :: r_ij, r_ik
      logical :: Zk1, Zk2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif

         n_neighbours_coordination = n_neighbours(at,i,max_dist=this%coordination_cutoff)

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij)
            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, n_neighbours(at,i)
               if( n == m ) cycle
               k = neighbour(at, i, m, distance = r_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)
               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 3 + n_neighbours_coordination
            enddo
         enddo
      enddo


   endsubroutine co_angle_3b_sizes

   subroutine co_distance_2b_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(co_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      real(dp) :: r_ij
      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at,i,n,distance=r_ij)
            if( r_ij >= this%cutoff ) cycle
!if( r_ij < 3.5_dp ) cycle


            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            n_descriptors = n_descriptors + 1
            n_cross = n_cross + 4 + n_neighbours(at,i,max_dist=this%coordination_cutoff) + n_neighbours(at,j,max_dist=this%coordination_cutoff)
         enddo
      enddo

   endsubroutine co_distance_2b_sizes

   subroutine cosnx_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(cosnx), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine cosnx_sizes

   subroutine trihis_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(trihis), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = at%N

      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_cross = n_cross + n_neighbours(at,i) + 1
      enddo

   endsubroutine trihis_sizes

   subroutine water_monomer_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(water_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(at%Z(i) == 8) then
            if(present(mask)) then
               if(.not. mask(i)) cycle
            endif
            n_descriptors = n_descriptors + 1
            n_cross = n_cross + 3
         endif
      enddo

   endsubroutine water_monomer_sizes

   subroutine water_dimer_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(water_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i, j, n
      real(dp) :: r_ij

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
call print("mask present ? "//present(mask)) 
      do i = 1, at%N
         if(at%Z(i) == 8) then
            if(present(mask)) then
               if(.not. mask(i)) cycle
            endif
            do n = 1, n_neighbours(at,i)
               j = neighbour(at,i,n,distance=r_ij)
               if(at%Z(j) == 8 .and. r_ij < this%cutoff) then
                  n_descriptors = n_descriptors + 1
                  n_cross = n_cross + 6
               endif
            enddo
         endif
      enddo
   endsubroutine water_dimer_sizes

   subroutine A2_dimer_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(A2_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i, j, iA1, iA2, iB1, iB2
      integer, dimension(at%N) :: A2_monomer_index
      real(dp) :: r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_sizes: descriptor object not initialised", error)
      endif

      call find_A2_monomer(at,this%atomic_number, this%monomer_cutoff, A2_monomer_index)

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         iA1 = i
         iA2 = neighbour(at,i,A2_monomer_index(i),distance=r_A1_A2)
         if( iA1 > iA2 ) cycle

         do j = i + 1, at%N
            iB1 = j
            iB2 = neighbour(at,j,A2_monomer_index(j),distance=r_B1_B2)
            if( iB1 > iB2 ) cycle

            r_A1_B1 = distance_min_image(at,iA1,iB1)
            r_A1_B2 = distance_min_image(at,iA1,iB2)

            r_A2_B1 = distance_min_image(at,iA2,iB1)
            r_A2_B2 = distance_min_image(at,iA2,iB2)
            
            if( all( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) < this%cutoff) ) then
               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 4
            endif
         enddo
      enddo

   endsubroutine A2_dimer_sizes

   subroutine AB_dimer_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(AB_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i, j, n_monomers, iA1, iA2, iB1, iB2
      integer, dimension(:,:), allocatable :: AB_monomer_index
      real(dp) :: r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_sizes: descriptor object not initialised", error)
      endif

      if( count(at%Z == this%atomic_number1) == count(at%Z == this%atomic_number2) ) then
         n_monomers = count(at%Z == this%atomic_number1)
      else
         RAISE_ERROR("AB_dimer_sizes: number of monomer atoms 1 ("//count(at%Z == this%atomic_number1)//") not equal to number of monomer atoms 2 ("//count(at%Z == this%atomic_number1)//")",error)
      endif

      allocate(AB_monomer_index(2,n_monomers))
      call find_AB_monomer(at,(/this%atomic_number1,this%atomic_number2/), this%monomer_cutoff, AB_monomer_index)

      n_descriptors = 0
      n_cross = 0

      do i = 1, n_monomers
         iA1 = AB_monomer_index(1,i)
         iB1 = AB_monomer_index(2,i)
         do j = i + 1, n_monomers
            iA2 = AB_monomer_index(1,j)
            iB2 = AB_monomer_index(2,j)

            r_A1_B1 = distance_min_image(at,iA1,iB1)
            r_A2_B2 = distance_min_image(at,iA2,iB2)

            r_A1_A2 = distance_min_image(at,iA1,iA2)
            r_B1_B2 = distance_min_image(at,iB1,iB2)

            r_A1_B2 = distance_min_image(at,iA1,iB2)
            r_A2_B1 = distance_min_image(at,iA2,iB1)
            
            if( all( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) < this%cutoff) ) then
               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 4
            endif
         enddo
      enddo

      deallocate(AB_monomer_index)

   endsubroutine AB_dimer_sizes

   subroutine bond_real_space_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(bond_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      type(atoms) :: at_copy
      integer :: i, j, k, n, m, shift_j(3)

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         n_descriptors = n_descriptors + n_neighbours(at, i, max_dist=this%bond_cutoff)

         do n = 1, n_neighbours(at, i)
            j = neighbour(at, i, n, shift=shift_j, max_dist=this%bond_cutoff)

            if(j == 0) cycle

            at_copy = at
            call add_atoms(at_copy, 0.5_dp * (at%pos(:,i) + at%pos(:,j) + matmul(at%lattice,shift_j)), 1)
            call calc_connect(at_copy)

            do m = 1, n_neighbours(at_copy, at%N + 1)
               k = neighbour(at_copy, at%N + 1, m, max_dist=this%cutoff)

               if(k == 0) cycle

               if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

               n_cross = n_cross + 1
            enddo

            call finalise(at_copy)
         enddo
      enddo

   endsubroutine bond_real_space_sizes

   subroutine atom_real_space_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(atom_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = at%N
      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff)*2 
      enddo

   endsubroutine atom_real_space_sizes

   subroutine power_so3_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(power_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine power_so3_sizes

   subroutine power_SO4_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(power_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine power_SO4_sizes

   subroutine soap_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(soap), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( .not. any( at%Z(i) == this%Z ) .and. .not. any(this%Z == 0) ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if(this%global) then
         n_descriptors = 1
      endif

   endsubroutine soap_sizes

   subroutine AN_monomer_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(AN_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AN_monomer: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + this%N
         if(.not.this%do_atomic) exit
      enddo

   endsubroutine AN_monomer_sizes

   subroutine general_monomer_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(general_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error
      integer, dimension(:,:), allocatable :: monomer_index
      integer :: i
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_monomer_sizes: descriptor object not initialised", error)
      endif

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_index,this%signature,associated_to_monomer,this%cutoff,this%atom_ordercheck,error)
      n_descriptors = size(monomer_index,2)
      n_cross=size(monomer_index)
      if(.not. all(associated_to_monomer)) then
         if(this%strict) then
            RAISE_ERROR("general_monomer_sizes: not all atoms assigned to a monomer", error)
         else
            call print("WARNING: general_monomer_sizes: not all atoms assigned to a monomer")
         endif
      endif

      deallocate(monomer_index)
      deallocate(associated_to_monomer)

   endsubroutine general_monomer_sizes

   subroutine com_dimer_sizes(this,at,n_descriptors,n_cross,mask,error)

      type(com_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index, monomer_pairs
      integer, dimension(:), allocatable :: pairs_diffs_map
      real(dp), dimension(:,:), allocatable :: mean_pos_diffs
      logical, dimension(:), allocatable :: associated_to_monomer
      logical :: double_count

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("com_dimer_sizes: descriptor object not initialised", error)
      endif

      double_count = .false.

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,this%signature_one,associated_to_monomer,this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%monomers_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
      else
        call find_general_monomer(at,monomer_two_index,this%signature_two,associated_to_monomer,this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         if(this%strict) then
            RAISE_ERROR("com_dimer_sizes: not all atoms assigned to a monomer", error)
         else
            call print("WARNING: com_dimer_sizes: not all atoms assigned to a monomer")
         endif
      endif

      if (this%mpifind) then
         call print("Using find_monomer_pairs_MPI", PRINT_NERD)
         call find_monomer_pairs_MPI(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,double_count,this%cutoff,error=error,use_com=.true.)
      else
         call find_monomer_pairs(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,double_count,this%cutoff,error=error,use_com=.true.)
      end if
      n_descriptors = size(pairs_diffs_map)
      n_cross=n_descriptors*(size(this%signature_one)+size(this%signature_two))
   
      deallocate(associated_to_monomer)
      deallocate(pairs_diffs_map)
      deallocate(monomer_pairs)
      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(mean_pos_diffs)

   endsubroutine com_dimer_sizes

   subroutine general_dimer_sizes(this,at,n_descriptors,n_cross,mask,error)

      type(general_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index, monomer_pairs
      integer, dimension(:), allocatable :: pairs_diffs_map
      real(dp), dimension(:,:), allocatable :: mean_pos_diffs
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_dimer_sizes: descriptor object not initialised", error)
      endif

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,this%signature_one,associated_to_monomer,this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%monomers_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
      else
        call find_general_monomer(at,monomer_two_index,this%signature_two,associated_to_monomer,this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         if(this%strict) then
            RAISE_ERROR("general_dimer_sizes: not all atoms assigned to a monomer", error)
         else
            call print("WARNING: general_dimer_sizes: not all atoms assigned to a monomer")
         endif
      endif

      if (this%mpifind) then
         call print("Using find_monomer_pairs_MPI", PRINT_NERD)
         call find_monomer_pairs_MPI(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,this%double_count,this%cutoff,error=error,use_com=this%use_com)
      else
         call find_monomer_pairs(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,this%monomers_identical,this%double_count,this%cutoff,error=error,use_com=this%use_com)
      end if
      n_descriptors = size(pairs_diffs_map)
      n_cross=n_descriptors*(size(this%signature_one)+size(this%signature_two))
   
      deallocate(associated_to_monomer)
      deallocate(pairs_diffs_map)
      deallocate(monomer_pairs)
      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(mean_pos_diffs)

   endsubroutine general_dimer_sizes

   subroutine general_trimer_sizes(this,at,n_descriptors,n_cross,mask,error)
 
      type(general_trimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error
      integer, dimension(:), allocatable ::  triplets_diffs_map
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index, monomer_three_index, monomer_triplets
      real(dp), dimension(:,:), allocatable :: triplets_diffs
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_trimer_sizes: descriptor object not initialised", error)
      endif

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,this%signature_one,associated_to_monomer,this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%one_two_identical) then
         allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
         monomer_two_index = monomer_one_index
      else
         call find_general_monomer(at,monomer_two_index,this%signature_two,associated_to_monomer,this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if
      if (this%one_three_identical) then
         allocate(monomer_three_index(size(monomer_one_index,1),size(monomer_one_index,2)))
         monomer_three_index = monomer_one_index
      else if (this%two_three_identical) then
         allocate(monomer_three_index(size(monomer_two_index,1),size(monomer_two_index,2)))
         monomer_three_index = monomer_two_index
      else
         call find_general_monomer(at,monomer_three_index,this%signature_three,associated_to_monomer,this%monomer_three_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         if(this%strict) then
            RAISE_ERROR("general_trimer_sizes: not all atoms assigned to a monomer", error)
         else
            call print("WARNING: general_trimer_sizes: not all atoms assigned to a monomer")
         endif
      endif
    
      if (this%use_com) then
         RAISE_ERROR("general_trimer_calc: use_com=T not implemented yet", error)
      end if
      if(this%mpifind) then
         call print("Using find_monomer_triplets_MPI", PRINT_NERD)
         call find_monomer_triplets_MPI(at,monomer_triplets,triplets_diffs,triplets_diffs_map,monomer_one_index,monomer_two_index,monomer_three_index,this%one_two_identical,this%one_three_identical,this%two_three_identical,this%cutoff,error,use_com=.false.)
      else
         call find_monomer_triplets(at,monomer_triplets,triplets_diffs,triplets_diffs_map,monomer_one_index,monomer_two_index,monomer_three_index,this%one_two_identical,this%one_three_identical,this%two_three_identical,this%cutoff,error)
      end if
      n_descriptors = size(triplets_diffs_map)
      n_cross=n_descriptors*(size(this%signature_one)+size(this%signature_two)+size(this%signature_three))

      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(monomer_three_index)
      if(allocated(monomer_triplets)) deallocate(monomer_triplets)
      if(allocated(triplets_diffs)) deallocate(triplets_diffs)
      deallocate(associated_to_monomer)

   endsubroutine general_trimer_sizes

   subroutine rdf_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(rdf), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("rdf_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine rdf_sizes

   subroutine as_distance_2b_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(as_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      real(dp) :: r_ij
      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("as_distance_2b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at,i,n,distance=r_ij)
            if( r_ij > this%max_cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            n_descriptors = n_descriptors + 1
            n_cross = n_cross + 4 + n_neighbours(at,i,max_dist=this%coordination_cutoff) + n_neighbours(at,j,max_dist=this%coordination_cutoff)
         enddo
      enddo

   endsubroutine as_distance_2b_sizes

   subroutine molecule_lo_d_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(molecule_lo_d), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error
      integer, dimension(:,:), allocatable :: monomer_index
      integer :: i
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("molecule_lo_d_sizes: descriptor object not initialised", error)
      endif

      n_descriptors=1
      n_cross = this%n_atoms

!!$      allocate(associated_to_monomer(at%N))
!!$      associated_to_monomer=.false.
!!$
!!$      call find_general_monomer(at,monomer_index,this%signature,associated_to_monomer,this%cutoff,this%atom_ordercheck,error)
!!$      n_descriptors = size(monomer_index,2)
!!$      n_cross=size(monomer_index)


   endsubroutine molecule_lo_d_sizes

   subroutine alex_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(alex), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("alex_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine alex_sizes

   subroutine distance_Nb_sizes(this,at,n_descriptors,n_cross,mask,error)
      type(distance_Nb), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, optional, intent(out) :: error

      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2
      real(dp) :: r_ij

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_Nb_sizes: descriptor object not initialised", error)
      endif

      call distance_Nb_calc_get_clusters(this,at,n_descriptors=n_descriptors,mask=mask,error=error)
      n_cross = n_descriptors * this%order

   endsubroutine distance_Nb_sizes

   function descriptor_n_permutations(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error

      integer :: descriptor_n_permutations, i

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_DISTANCE_2b,DT_COORDINATION, &
            DT_ANGLE_3B,DT_CO_ANGLE_3B,DT_CO_DISTANCE_2b,DT_COSNX,DT_TRIHIS,DT_WATER_MONOMER,DT_BOND_REAL_SPACE,DT_ATOM_REAL_SPACE,DT_POWER_SO3,DT_POWER_SO4,DT_SOAP,DT_RDF, DT_ALEX, DT_COM_DIMER)

            descriptor_n_permutations = 1
            
         case(DT_WATER_DIMER)
            descriptor_n_permutations = NP_WATER_DIMER
         case(DT_A2_DIMER)
            descriptor_n_permutations = NP_A2_DIMER
         case(DT_AB_DIMER)
            descriptor_n_permutations = NP_AB_DIMER
         case(DT_AN_MONOMER)
            if(this%descriptor_AN_monomer%do_atomic) then
               descriptor_n_permutations = factorial(this%descriptor_AN_monomer%N-1)
            else
               descriptor_n_permutations = factorial(this%descriptor_AN_monomer%N)
            endif
         case(DT_GENERAL_MONOMER)
            if (.not. this%descriptor_general_monomer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_n_permutations: permutation_data not initialised "//this%descriptor_type,error)
            end if
            descriptor_n_permutations = this%descriptor_general_monomer%permutation_data%n_perms
         case(DT_GENERAL_DIMER)
            if (.not. this%descriptor_general_dimer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_n_permutations: permutation_data not initialised "//this%descriptor_type,error)
            end if
            descriptor_n_permutations = this%descriptor_general_dimer%permutation_data%n_perms
         case(DT_GENERAL_TRIMER)
            if (.not. this%descriptor_general_trimer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_n_permutations: permutation_data not initialised "//this%descriptor_type,error)
            end if
            descriptor_n_permutations = this%descriptor_general_trimer%permutation_data%n_perms
         case(DT_MOLECULE_LO_D)
            if (.not. this%descriptor_molecule_lo_d%permutation_data%initialised)then
              RAISE_ERROR("descriptor_n_permutations: permutation_data not initialised "//this%descriptor_type,error)
            end if
            descriptor_n_permutations = this%descriptor_molecule_lo_d%permutation_data%n_perms
         case(DT_DISTANCE_NB)
            descriptor_n_permutations = this%descriptor_distance_Nb%n_permutations
         case default
            RAISE_ERROR("descriptor_n_permutations: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_n_permutations

   subroutine descriptor_permutations(this,permutations,error)
      type(descriptor), intent(in) :: this
      type(permutation_data_type) :: my_permutation_data
      integer, dimension(:,:), intent(out) :: permutations
      integer, optional, intent(out) :: error

      integer :: i, d, np, n, m, ip, j
      integer,dimension(1) :: unit_vec
      integer, dimension(:), allocatable :: this_perm
      integer, dimension(:,:), allocatable :: distance_matrix, atom_permutations, sliced_permutations

      INIT_ERROR(error)

      d = descriptor_dimensions(this,error)
      np = descriptor_n_permutations(this,error)
      call check_size('permutations',permutations, (/d,np/),'descriptor_permutations',error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_DISTANCE_2b,DT_COORDINATION, &
            DT_ANGLE_3B,DT_CO_ANGLE_3B,DT_CO_DISTANCE_2b,DT_COSNX,DT_TRIHIS,DT_WATER_MONOMER,DT_BOND_REAL_SPACE,DT_ATOM_REAL_SPACE,DT_POWER_SO3,DT_POWER_SO4,DT_SOAP,DT_RDF, DT_ALEX, DT_COM_DIMER)
            
            permutations(:,1) = (/ (i, i = 1, size(permutations,1)) /)
         case(DT_WATER_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15/) ! original order
            permutations(:,2) = (/1, 3, 2, 4, 5, 7, 6, 8, 9, 10, 13, 14, 11, 12, 15/) ! swap Hs on monomer A
            permutations(:,3) = (/1, 2, 3, 5, 4, 6, 7, 9, 8, 10, 12, 11, 14, 13, 15/) ! swap Hs on monomer B
            permutations(:,4) = (/1, 3, 2, 5, 4, 7, 6, 9, 8, 10, 14, 13, 12, 11, 15/) ! swap Hs on both monomers
            permutations(:,5) = (/1, 8, 9, 6, 7, 4, 5, 2, 3, 15, 11, 13, 12, 14, 10/) ! swap monomers A and B
            permutations(:,6) = (/1, 9, 8, 6, 7, 5, 4, 2, 3, 15, 12, 14, 11, 13, 10/) ! swap monomers and Hs on monomer A
            permutations(:,7) = (/1, 8, 9, 7, 6, 4, 5, 3, 2, 15, 13, 11, 14, 12, 10/) ! swap monomers and Hs on monomer B
            permutations(:,8) = (/1, 9, 8, 7, 6, 5, 4, 3, 2, 15, 14, 12, 13, 11, 10/) ! swap monomers and Hs on both monomers

         case(DT_A2_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6/) ! original order
            permutations(:,2) = (/1, 2, 5, 6, 3, 4/) ! swap atoms on monomer A
            permutations(:,3) = (/1, 2, 4, 3, 6, 5/) ! swap atoms on monomer B
            permutations(:,4) = (/1, 2, 6, 5, 4, 3/) ! swap atoms on both monomers
            permutations(:,5) = (/2, 1, 3, 5, 4, 6/) ! swap monomers A and B
            permutations(:,6) = (/2, 1, 5, 3, 6, 4/) ! swap monomers and atoms on monomer A
            permutations(:,7) = (/2, 1, 4, 6, 3, 5/) ! swap monomers and atoms on monomer B
            permutations(:,8) = (/2, 1, 6, 4, 5, 3/) ! swap monomers and atoms on both monomers
            
         case(DT_AB_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6/) ! original order
            permutations(:,2) = (/2, 1, 3, 4, 6, 5/) ! swap monomers

         case(DT_AN_MONOMER)
            allocate(distance_matrix(this%descriptor_AN_monomer%N,this%descriptor_AN_monomer%N), atom_permutations(this%descriptor_AN_monomer%N,np))

            if(this%descriptor_AN_monomer%do_atomic) then
               atom_permutations(1,:) = 0
               call generate_AN_permutations(atom_permutations(2:this%descriptor_AN_monomer%N,:))
               atom_permutations = atom_permutations + 1
            else
               call generate_AN_permutations(atom_permutations(:,:))
            endif

            i = 0
            distance_matrix = 0
            do n = 2, this%descriptor_AN_monomer%N
               i = i + 1
               distance_matrix(1,n) = i
               distance_matrix(n,1) = i
               do m = n+1, this%descriptor_AN_monomer%N
                  i = i + 1
                  distance_matrix(m,n) = i
                  distance_matrix(n,m) = i
               enddo
            enddo

            do ip = 1, np
               i = 0
               do n = 2, this%descriptor_AN_monomer%N
                  i = i + 1
                  permutations(i,ip) = distance_matrix(atom_permutations(1,ip),atom_permutations(n,ip))
                  do m = n+1, this%descriptor_AN_monomer%N
                     i = i + 1
                     permutations(i,ip) = distance_matrix(atom_permutations(m,ip),atom_permutations(n,ip))
                  enddo
               enddo
            enddo
            deallocate(distance_matrix,atom_permutations)


         case(DT_GENERAL_MONOMER)
            if (.not. this%descriptor_general_monomer%permutation_data%initialised) then
              RAISE_ERROR("descriptor_permutations: permutation_data not initialised "//this%descriptor_type,error)
            else if (this%descriptor_general_monomer%permutation_data%perm_number /= 1) then
              RAISE_ERROR("descriptor_permutations: permutation_data%perm_number must be initialised to one"//this%descriptor_type,error)
            end if

            call permutation_data_copy(my_permutation_data, this%descriptor_general_monomer%permutation_data)

            if (my_permutation_data%n_perms > 1) then
              call next(my_permutation_data, 1)
            end if

            permutations=my_permutation_data%dist_vec_permutations

         case(DT_GENERAL_DIMER)
            if (.not. this%descriptor_general_dimer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_permutations: permutation_data not initialised "//this%descriptor_type,error)
            else if (this%descriptor_general_dimer%permutation_data%perm_number /= 1) then
              RAISE_ERROR("descriptor_permutations: permutation_data%perm_number must be initialised to one"//this%descriptor_type,error)
            end if

            call permutation_data_copy(my_permutation_data, this%descriptor_general_dimer%permutation_data)

            if (my_permutation_data%n_perms > 1) then
              call next(my_permutation_data, 1)
            end if

            permutations=my_permutation_data%dist_vec_permutations

         case(DT_GENERAL_TRIMER)
            if (.not. this%descriptor_general_trimer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_permutations: permutation_data not initialised "//this%descriptor_type,error)
            else if (this%descriptor_general_trimer%permutation_data%perm_number /= 1) then
              RAISE_ERROR("descriptor_permutations: permutation_data%perm_number must be initialised to one"//this%descriptor_type,error)
            end if

            call permutation_data_copy(my_permutation_data, this%descriptor_general_trimer%permutation_data)

            if (my_permutation_data%n_perms > 1) then
              call next(my_permutation_data, 1)
            end if

            permutations=my_permutation_data%dist_vec_permutations

         case(DT_MOLECULE_LO_D)
            if (.not. this%descriptor_molecule_lo_d%permutation_data%initialised) then
              RAISE_ERROR("descriptor_permutations: permutation_data not initialised "//this%descriptor_type,error)
            else if (this%descriptor_molecule_lo_d%permutation_data%perm_number /= 1) then
              RAISE_ERROR("descriptor_permutations: permutation_data%perm_number must be initialised to one"//this%descriptor_type,error)
            end if

            call permutation_data_copy(my_permutation_data, this%descriptor_molecule_lo_d%permutation_data)

            if (my_permutation_data%n_perms > 1) then
              call next(my_permutation_data, 1)
            end if

            allocate(sliced_permutations(size(this%descriptor_molecule_lo_d%included_components),my_permutation_data%n_perms))
            allocate(this_perm(size(this%descriptor_molecule_lo_d%included_components)))
            sliced_permutations =my_permutation_data%dist_vec_permutations(this%descriptor_molecule_lo_d%included_components,:)

            do j=1,my_permutation_data%n_perms
              this_perm=sliced_permutations(:,j)
              do i=1,size(this%descriptor_molecule_lo_d%included_components)
                unit_vec=maxloc(this%descriptor_molecule_lo_d%included_components, mask=this%descriptor_molecule_lo_d%included_components .eq. this_perm(i))
                if (unit_vec(1) == 0) then
                  RAISE_ERROR("descriptor_permutations: you have specified symmetries between atoms with different connectivity",error)
                end if
                permutations(i,j) =unit_vec(1)
              end do
            end do
! begin brau
            if(size(this%descriptor_molecule_lo_d%included_components) > maxval(this%descriptor_molecule_lo_d%included_components)) then
              permutations=my_permutation_data%dist_vec_permutations
            end if
! end brau
         case(DT_DISTANCE_NB)
            permutations = this%descriptor_distance_Nb%permutations
         case default
            RAISE_ERROR("descriptor_permutations: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_permutations

   subroutine generate_AN_permutations(this,list,error)
      integer, dimension(:,:), intent(out) :: this
      integer, dimension(:), intent(in), optional :: list
      integer, optional, intent(out) :: error

      integer, dimension(:), allocatable :: my_list, my_list_uniq
      integer :: i, n, m, p, np, min_tail, min_tail_i, tmp_i

      INIT_ERROR(error)

      if(present(list)) then
         n = size(list)
         allocate(my_list(n))
         my_list = list
      else
         n = size(this,1)
         allocate(my_list(n))
         my_list = (/(i, i = 1, n)/)
      endif

      call uniq(my_list, my_list_uniq)

      np = factorial(size(my_list_uniq))

      call check_size('this', this, (/n,np/), 'generate_permutations',error)

      call sort_array(my_list)

      this(:,1) = my_list

      do p = 2, np
         ! Find longest tail that is ordered in decreasing order.
         do m = n - 1, 1, -1
            if(my_list(m) < my_list(m+1)) exit
         enddo

         min_tail = my_list(m+1)
         min_tail_i = m+1
         ! Find the smallest number bigger than my_list(m) in the tail
         do i = m + 1, n
            if(min_tail > my_list(i) .and. my_list(m) < my_list(i)) then
               min_tail = my_list(i)
               min_tail_i = i
            endif
         enddo

         ! swap
         tmp_i = my_list(m)
         my_list(m) = my_list(min_tail_i)
         my_list(min_tail_i) = tmp_i

         
         ! reverse tail
         my_list(m+1:n) = my_list(n:m+1:-1)

         this(:,p) = my_list
      enddo

   endsubroutine generate_AN_permutations

   subroutine real_space_fourier_coefficients(at,l_max,atom_coefficient)
      type(atoms), intent(in) :: at
      integer, intent(in) :: l_max
      type(neighbour_type), dimension(:), allocatable :: atom_coefficient

      integer :: i, j, n, l, m
      real(dp) :: r
      real(dp), dimension(3) :: d

      if(.not.allocated(atom_coefficient)) allocate(atom_coefficient(at%N))

      do i = 1, at%N
         if(.not. allocated(atom_coefficient(i)%neighbour)) allocate(atom_coefficient(i)%neighbour(n_neighbours(at,i)))
         do n = 1, n_neighbours(at,i)

            j = neighbour(at,i,n,distance = r, diff = d)
            atom_coefficient(i)%neighbour(n)%r = r
            atom_coefficient(i)%neighbour(n)%u = d / r

            if(.not. allocated(atom_coefficient(i)%neighbour(n)%spherical_harmonics)) allocate( atom_coefficient(i)%neighbour(n)%spherical_harmonics(0:l_max), &
            atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(0:l_max) )
            do l = 0, l_max
               if(.not. allocated(atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m)) &
               allocate(atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m(-l:l))
               if(.not. allocated(atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm)) &
               allocate(atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm(3,-l:l))

               atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m = CPLX_ZERO
               atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm = CPLX_ZERO

               do m = -l, l
                  atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m(m) = SphericalYCartesian(l,m,d)
                  atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm(:,m) = GradSphericalYCartesian(l,m,d)
               enddo
            enddo
         enddo
      enddo

   endsubroutine real_space_fourier_coefficients

   function real_space_covariance_coefficient(anc1,anc2,i1,i2,alpha,l_max,f1,f2)
      type(neighbour_type), dimension(:), intent(in) :: anc1, anc2
      real(dp), intent(in) :: alpha
      integer, intent(in) :: i1, i2, l_max
      real(dp), dimension(:,:), intent(out), optional :: f1, f2

      real(dp) :: real_space_covariance_coefficient

      complex(dp) :: real_space_covariance_in, I_lm1m2
      integer :: n1, n2, l, m1, m2, k
      real(dp) :: r1, r2, arg_bess, fac_exp, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_lp, grad_mo_spher_bess_fi_ki_l
      real(dp), dimension(3) :: u1, u2, grad_arg_bess1, grad_fac_exp1, grad_arg_bess2, grad_fac_exp2
      type(cplx_2d), dimension(:), allocatable :: integral_r
      type(grad_spherical_harmonics_overlap_type), dimension(:), allocatable :: grad_integral_r1, grad_integral_r2

      logical :: do_derivative

      do_derivative = (present(f1) .or. present(f2)) 

      real_space_covariance_in = CPLX_ZERO

      allocate(integral_r(0:l_max))
      do l = 0, l_max
         allocate(integral_r(l)%mm(-l:l,-l:l))
         integral_r(l)%mm = CPLX_ZERO
      enddo

      if(present(f1)) then
         allocate(grad_integral_r1(0:size(anc1(i1)%neighbour)))
         do n1 = 0, size(anc1(i1)%neighbour)
            allocate(grad_integral_r1(n1)%grad_integral(0:l_max))
            do l = 0, l_max
               allocate(grad_integral_r1(n1)%grad_integral(l)%mm(3,-l:l,-l:l))
               grad_integral_r1(n1)%grad_integral(l)%mm = CPLX_ZERO
            enddo
         enddo
      endif

      if(present(f2)) then
         allocate(grad_integral_r2(0:size(anc2(i2)%neighbour)))
         do n2 = 0, size(anc2(i2)%neighbour)
            allocate(grad_integral_r2(n2)%grad_integral(0:l_max))
            do l = 0, l_max
               allocate(grad_integral_r2(n2)%grad_integral(l)%mm(3,-l:l,-l:l))
               grad_integral_r2(n2)%grad_integral(l)%mm = CPLX_ZERO
            enddo
         enddo
      endif
      do n1 = 1, size(anc1(i1)%neighbour)
         r1 = anc1(i1)%neighbour(n1)%r
         u1 = anc1(i1)%neighbour(n1)%u
         do n2 = 1, size(anc2(i2)%neighbour)
            r2 = anc2(i2)%neighbour(n2)%r

            u2 = anc2(i2)%neighbour(n2)%u

            arg_bess = alpha*r1*r2
            fac_exp = exp(-0.5_dp*alpha*(r1**2+r2**2))

            if(present(f1)) then
               grad_arg_bess1 = alpha*r2*u1
               grad_fac_exp1 = -fac_exp*alpha*r1*u1
            endif

            if(present(f2)) then
               grad_arg_bess2 = alpha*r1*u2
               grad_fac_exp2 = -fac_exp*alpha*r2*u2
            endif

            do l = 0, l_max
               if( l == 0 ) then
                  mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                  mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                  if(do_derivative) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
               else
                  mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                  mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                  if(do_derivative) then
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                     mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                  else
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess
                  endif

               endif


               if(do_derivative) grad_mo_spher_bess_fi_ki_l = 0.5_dp * (mo_spher_bess_fi_ki_lp - mo_spher_bess_fi_ki_l / arg_bess + mo_spher_bess_fi_ki_lm)
                  
               do m1 = -l, l
                  do m2 = -l, l
                     I_lm1m2 = conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * mo_spher_bess_fi_ki_l*fac_exp
                     integral_r(l)%mm(m2,m1) = integral_r(l)%mm(m2,m1) + I_lm1m2
                     if(present(f1)) then
                        grad_integral_r1(n1)%grad_integral(l)%mm(:,m2,m1) = grad_integral_r1(n1)%grad_integral(l)%mm(:,m2,m1) + &
                        anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * &
                        ( conjg(anc1(i1)%neighbour(n1)%grad_spherical_harmonics(l)%mm(:,m1)) * mo_spher_bess_fi_ki_l*fac_exp + &
                        conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * ( grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp1 ) )
                     endif

                     if(present(f2)) then
                        grad_integral_r2(n2)%grad_integral(l)%mm(:,m2,m1) = grad_integral_r2(n2)%grad_integral(l)%mm(:,m2,m1) + &
                        conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * &
                        ( anc2(i2)%neighbour(n2)%grad_spherical_harmonics(l)%mm(:,m2) * mo_spher_bess_fi_ki_l*fac_exp + &
                        anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * ( grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp2 ) )
                     endif

                  enddo
               enddo
            enddo
         enddo
      enddo

      if(present(f1)) then
         f1 = 0.0_dp
         do n1 = 0, size(anc1(i1)%neighbour)
            do l = 0, l_max
               do k = 1, 3
                  f1(k,n1+1) = f1(k,n1+1) + real(sum(conjg(grad_integral_r1(n1)%grad_integral(l)%mm(k,:,:))*integral_r(l)%mm(:,:)))
               enddo
            enddo
         enddo
         f1 = 2.0_dp * f1
      endif

      if(present(f2)) then
         f2 = 0.0_dp
         do n2 = 0, size(anc2(i2)%neighbour)
            do l = 0, l_max
               do k = 1, 3
                  f2(k,n2+1) = f2(k,n2+1) + real(sum(conjg(grad_integral_r2(n2)%grad_integral(l)%mm(k,:,:))*integral_r(l)%mm(:,:)))
               enddo
            enddo
         enddo
         f2 = 2.0_dp * f2
      endif

      do l = 0, l_max
         real_space_covariance_in = real_space_covariance_in + sum(conjg(integral_r(l)%mm) * integral_r(l)%mm)
      enddo
      real_space_covariance_coefficient = real(real_space_covariance_in)

      do l = 0, l_max
         deallocate(integral_r(l)%mm)
      enddo
      deallocate(integral_r)

      if(present(f1)) then
         do n1 = 0, size(anc1(i1)%neighbour)
            do l = 0, l_max
               deallocate(grad_integral_r1(n1)%grad_integral(l)%mm)
            enddo
            deallocate(grad_integral_r1(n1)%grad_integral)
         enddo
         deallocate(grad_integral_r1)
      endif

      if(present(f2)) then
         do n2 = 0, size(anc2(i2)%neighbour)
            do l = 0, l_max
               deallocate(grad_integral_r2(n2)%grad_integral(l)%mm)
            enddo
            deallocate(grad_integral_r2(n2)%grad_integral)
         enddo
         deallocate(grad_integral_r2)
      endif

   endfunction real_space_covariance_coefficient

   function real_space_covariance(at1,at2,i1,i2,alpha,l_max,f1,f2)
      type(atoms), intent(in) :: at1, at2
      real(dp), intent(in) :: alpha
      integer, intent(in) :: i1, i2, l_max
      real(dp), dimension(:,:), intent(inout), optional :: f1, f2

      real(dp) :: real_space_covariance

      complex(dp) :: real_space_covariance_in, I_lm1m2
      integer :: j1, j2, n1, n2, l, m1, m2
      real(dp) :: r1, r2, arg_bess, fac_exp, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm
      real(dp), dimension(3) :: d1, d2
      type(cplx_2d), dimension(:), allocatable :: integral_r

      logical :: do_derivative

      do_derivative = (present(f1) .or. present(f2)) 

      real_space_covariance_in = CPLX_ZERO

      allocate(integral_r(0:l_max))
      do l = 0, l_max
         allocate(integral_r(l)%mm(-l:l,-l:l))
         integral_r(l)%mm = CPLX_ZERO
      enddo

      do n1 = 1, n_neighbours(at1,i1)
         j1 = neighbour(at1,i1,n1,distance = r1, diff = d1)
         do n2 = 1, n_neighbours(at2,i2)
            j2 = neighbour(at2,i2,n2,distance = r2, diff = d2)

            arg_bess = alpha*r1*r2
            fac_exp = exp(-0.5_dp*alpha*(r1**2+r2**2))

            do l = 0, l_max
               if( l == 0 ) then
                  mo_spher_bess_fi_ki_lmm = sinh(arg_bess)/arg_bess
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm
               elseif( l == 1 ) then
                  mo_spher_bess_fi_ki_lm = ( arg_bess*cosh(arg_bess) - sinh(arg_bess) ) / arg_bess**2
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lm
               else
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l+1)*mo_spher_bess_fi_ki_lm / arg_bess
                  mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                  mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
               endif
                  
               do m1 = -l, l
                  do m2 = -l, l
                     I_lm1m2 = conjg(SphericalYCartesian(l,m1,d1)) * SphericalYCartesian(l,m2,d2)*mo_spher_bess_fi_ki_l*fac_exp
                     integral_r(l)%mm(m2,m1) = integral_r(l)%mm(m2,m1) + I_lm1m2
                  enddo
               enddo
            enddo
         enddo
      enddo

      do l = 0, l_max
         real_space_covariance_in = real_space_covariance_in + sum(conjg(integral_r(l)%mm) * integral_r(l)%mm)
      enddo
      real_space_covariance = real(real_space_covariance_in)

      do l = 0, l_max
         deallocate(integral_r(l)%mm)
      enddo
      deallocate(integral_r)

   endfunction real_space_covariance

   function RadialFunction(this,r,i)
      type(RadialFunction_type), intent(in) :: this
      real(dp), intent(in) :: r
      integer, intent(in) :: i

      real(dp) :: RadialFunction
   
      real(dp), dimension(this%n_max) :: h
      integer :: j
   
      if( r < this%cutoff ) then
         do j = 1, this%n_max
            h(j) = (this%cutoff-r)**(j+2) / this%NormFunction(j)
         enddo
         RadialFunction = dot_product(this%RadialTransform(:,i),h)
      else
         RadialFunction = 0.0_dp
      endif
   
   endfunction RadialFunction
   
   function GradRadialFunction(this,r,i)
      type(RadialFunction_type), intent(in) :: this
      real(dp), intent(in) :: r
      integer, intent(in) :: i

      real(dp) :: GradRadialFunction
   
      real(dp), dimension(this%n_max) :: h
      integer :: j
   
      if( r < this%cutoff ) then
         do j = 1, this%n_max
            h(j) = - (j+2) * (this%cutoff-r)**(j+1) / this%NormFunction(j)
         enddo
         GradRadialFunction = dot_product(this%RadialTransform(:,i),h)
      else
         GradRadialFunction = 0.0_dp
      endif
   
   endfunction GradRadialFunction

   subroutine bond_list_next_layer(shallow_list,deep_list,error)
      ! shallow list contains all pairs of bonded atoms, whereas deep_list will include pairs of atoms
      ! separated by two bonds, three bonds, etc., depending on how many times this function has been called
      type(Table), intent(in):: shallow_list
      type(Table), intent(inout):: deep_list
      type(Table) :: connected, pairs_containing_atom_i, pairs_containing_atom_j, deep_list_input, bonded_to_atom_i, bonded_to_atom_j
      integer :: i, atom_i, atom_j, j, atom_k, k, N_input,N_deep,N_shallow
      integer, intent(inout), optional :: error
      logical, dimension(:), allocatable :: mask_deep, mask_shallow
      logical :: i_k_present, k_i_present, j_k_present, k_j_present

      allocate(mask_deep(deep_list%N))
      allocate(mask_shallow(shallow_list%N))
      mask_deep = .False.
      mask_shallow = .False.

      ! make a copy of the deep list input
      deep_list_input = deep_list
      N_input=deep_list_input%N
      N_shallow=shallow_list%N

      ! loop over pairs in deep_list_input
      do i=1,N_input


        atom_i = deep_list_input%int(1,i)
        atom_j = deep_list_input%int(2,i)

        ! select the 1st neighbours of atom_i and atom_j
        mask_shallow = shallow_list%int(1,:N_shallow) .eq. atom_i .or. shallow_list%int(2,:N_shallow) .eq. atom_i
        call select(bonded_to_atom_i, shallow_list, row_mask=mask_shallow)

        mask_shallow = shallow_list%int(1,:N_shallow) .eq. atom_j .or. shallow_list%int(2,:N_shallow) .eq. atom_j
        call select(bonded_to_atom_j, shallow_list, row_mask=mask_shallow)


        do k=1,bonded_to_atom_j%N ! and append distances between atom i and all atoms bonded to atom j
          N_deep = deep_list%N
          atom_k = bonded_to_atom_j%int(1,k)
          if (atom_k .eq. atom_j) then
            atom_k = bonded_to_atom_j%int(2,k) 
          end if
          if (atom_k .eq. atom_i) cycle

          i_k_present = any(deep_list%int(1,:N_deep) .eq. atom_i .and. (deep_list%int(2,:N_deep) .eq. atom_k))
          k_i_present = any(deep_list%int(1,:N_deep) .eq. atom_k .and. (deep_list%int(2,:N_deep) .eq. atom_i))

          if (.not. k_i_present .and. .not. i_k_present) then
            call append(deep_list,(/atom_i,atom_k/))
          end if
        end do

        do k=1,bonded_to_atom_i%N ! and append distances between atom i and all atoms bonded to atom j
          N_deep = deep_list%N
          atom_k = bonded_to_atom_i%int(1,k)
          if (atom_k .eq. atom_i) then
            atom_k = bonded_to_atom_i%int(2,k) 
          end if
          if (atom_k .eq. atom_j) cycle

          j_k_present = any(deep_list%int(1,:N_deep) .eq. atom_j .and. (deep_list%int(2,:N_deep) .eq. atom_k))
          k_j_present = any(deep_list%int(1,:N_deep) .eq. atom_k .and. (deep_list%int(2,:N_deep) .eq. atom_j))

          if (.not. k_j_present .and. .not. j_k_present) then
            call append(deep_list,(/atom_j,atom_k/))
          end if
        end do

      end do


   end subroutine bond_list_next_layer

   subroutine transfer_initialise(this, args_str, error)
      type(transfer_parameters_type), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)
      call initialise(params)
      call param_register(params, 'do_transfer', 'false', this%do_transfer, help_string="Enable transfer function")
      call param_register(params, 'transfer_factor', '5.0', this%factor, help_string="Transfer function: stretch factor")
      call param_register(params, 'transfer_width', '1.0', this%width, help_string="Transfer function: transition width")
      call param_register(params, 'transfer_r0', '3.0', this%r0, help_string="Transfer function: transition distance")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true., task='transfer_initialise args_str')) then
         RAISE_ERROR("transfer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if (this%do_transfer) then
         call print("Using transfer function with factor="//this%factor//", r0="//this%r0//", width="//this%width)
      endif
   end subroutine transfer_initialise

   function transferfunction_grad(x, params) result(td)
     real(dp), intent(in) :: x
     type(transfer_parameters_type), intent(in) :: params
     real(dp) :: td
     ! for x << r0 - width: params%factor
     ! for x >> r0 + width: 1
     td = (params%factor - 1.0_dp) * 0.5_dp*(tanh((params%r0 - x)/params%width) + 1.0_dp) + 1.0_dp
   end function transferfunction_grad

   function transferfunction(x, params) result(t)
     real(dp), intent(in) :: x
     type(transfer_parameters_type), intent(in) :: params
     real(dp) :: t
     ! for x >> r0 - width: identity (slope=1)
     ! for x << r0 + width: linear slope = factor
     t = (1.0_dp - params%factor) * 0.5_dp*(params%width*(log(2.0_dp*cosh((x - params%r0)/params%width))) + x + params%r0) + params%factor * x
   end function transferfunction

   function graphIsConnected(connectivityMatrix,error)

      logical, dimension(:,:), intent(in) :: connectivityMatrix
      integer, intent(out), optional :: error
      logical :: graphIsConnected

      logical, dimension(:), allocatable :: visitedVertices

      INIT_ERROR(error)

      if( .not. is_square(connectivityMatrix) ) then
         RAISE_ERROR("graphIsConnected: not square matrix",error)
      endif

      allocate(visitedVertices(size(connectivityMatrix,1)))

      call graphBFS(connectivityMatrix,1,visitedVertices=visitedVertices,error=error)
      graphIsConnected = all(visitedVertices)

      deallocate(visitedVertices)

   endfunction graphIsConnected

   subroutine graphBFS(connectivityMatrix,startVertex,visitedVertices,tree,error)

      logical, dimension(:,:), intent(in) :: connectivityMatrix
      integer, intent(in) :: startVertex
      logical, dimension(:), target, intent(out), optional :: visitedVertices
      integer, dimension(:,:), allocatable, intent(out), optional :: tree
      integer, intent(out), optional :: error

      type(LinkedList_i1d), pointer :: LL_edges => null(), LL_remove => null(), LL_tree => null()

      logical, dimension(:), pointer :: my_visitedVertices
      integer, dimension(:), pointer :: edge
      integer, dimension(2) :: vw

      INIT_ERROR(error)

      if( .not. is_square(connectivityMatrix) ) then
         RAISE_ERROR("graphBFS: not square matrix",error)
      endif

      if( present( visitedVertices ) ) then
         my_visitedVertices => visitedVertices
      else
         allocate(my_visitedVertices(size(connectivityMatrix,1)))
      endif

      my_visitedVertices = .false.
      call graphSearch(connectivityMatrix,startVertex,LL_edges,my_visitedVertices,error)
      do while( associated(LL_edges) )
         LL_remove => LL_edges
         edge => retrieve_node(LL_remove)
         vw = edge
         call delete_node(LL_edges,LL_remove)
         if( .not. my_visitedVertices(vw(2)) ) then

            if(present(tree)) call append(LL_tree,vw)
            call graphSearch(connectivityMatrix, vw(2), LL_edges, my_visitedVertices,error)
         endif
      enddo

      if( .not. present( visitedVertices ) ) deallocate(my_visitedVertices)

      if (present(tree)) then
         call retrieve(LL_tree,tree)
         call finalise(LL_tree)
      endif

   endsubroutine graphBFS

   subroutine graphSearch(connectivityMatrix, vertex, LL_edges, visitedVertices,error)
      logical, dimension(:,:), intent(in) :: connectivityMatrix
      integer, intent(in) :: vertex
      type(LinkedList_i1d), pointer, intent(inout) :: LL_edges
      logical, dimension(:), intent(inout) :: visitedVertices
      integer, intent(out), optional :: error

      integer :: i

      INIT_ERROR(error)

      if( .not. is_square(connectivityMatrix) ) then
         RAISE_ERROR("graphSearch: not square matrix",error)
      endif

      visitedVertices(vertex) = .true.

      do i = 1, size(connectivityMatrix,1)
         if( connectivityMatrix(i,vertex) ) call append(LL_edges,(/vertex,i/))
      enddo

   endsubroutine graphSearch

endmodule descriptors_module
