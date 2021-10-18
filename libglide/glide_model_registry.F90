!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   glide_model_registry.F90 - part of the Community Ice Sheet Model (CISM)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2018
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_model_registry

  use glimmer_log
  use glide_types

  implicit none
  private

  !> module for counting and registering all models (i.e., model instances for different
  !> ice sheets) included in this run

  ! Usage is as follows:
  !
  ! Required:
  ! - Call register_model for each model in initialization (this should be called
  !   relatively early in initialization, so that model%model_id is set before it's
  !   needed)
  ! - Call deregister_model for each model in finalization
  !
  ! Optional:
  ! - Call set_max_models once in initialization. If this is called, it should be called
  !   before any calls to register_model. It is optional to call this; if it is not
  !   called, things will still work correctly, but some arrays will be allocated to be
  !   larger than they need to be (based on the value of default_max_models).
  ! - Call check_num_models once in initialization. It only makes sense to call this if
  !   set_max_models was called; if it is called, it should be called after all calls to
  !   register_model. This checks that the number of calls to register_model matched the
  !   expected number as set in the call to set_max_models.

  ! default_max_models is used if set_max_models is never called
  integer, parameter :: default_max_models = 8
  integer :: max_models = default_max_models
  integer :: num_models = 0  ! number of models actually registered so far

  type pmodel_type
     !> Contains a pointer to a model
     !> This is a hack to get around Fortran's lack of arrays of pointers
     type(glide_global_type), pointer :: p => null()
  end type pmodel_type

  !> Pointers to all registered models
  type(pmodel_type), allocatable, public, save :: registered_models(:)

  public :: set_max_models
  public :: register_model
  public :: check_num_models
  public :: deregister_model
  public :: get_num_models
  public :: get_max_models

contains

  subroutine set_max_models(max_models_in)
    !> Set the maximum number of models that will be used in this run.
    !>
    !> If this is called, it should be called before any calls to register_model. It is
    !> optional to call this; if it is not called, things will still work correctly, but
    !> some arrays will be allocated to be larger than they need to be (based on the value
    !> of default_max_models).
    integer, intent(in) :: max_models_in

    if (num_models > 0) then
       call write_log("Attempt to call set_max_models after register_model has already been called", &
            GM_FATAL)
    end if

    max_models = max_models_in
  end subroutine set_max_models

  subroutine register_model(model)
    !> Registers a model, setting model%model_id and ensuring that it is finalised in the
    !> case of an error.
    !>
    !> This should be called relatively early in initialization so that model%model_id is
    !> set before it's needed.
    type(glide_global_type), target :: model

    num_models = num_models + 1

    ! num_models == 1 indicates that this is the first time register_model is being
    ! called. This implies that registered_models is not yet allocated, so we need to
    ! allocate it.
    if (num_models == 1) then
       if (allocated(registered_models)) then
          ! This most likely indicates a programming error, not a user error
          call write_log("Somehow registered_models is already allocated", GM_FATAL)
       end if

       allocate(registered_models(max_models))
    end if

    if (num_models > max_models) then
       ! This can happen if set_max_models gave the wrong number, or (probably more
       ! likely) if set_max_models wasn't called and so we are using the default value
       ! for max_models; in the latter case, the code should either be changed to call
       ! set_max_models, or the hard-coded default_max_models should be increased.
       call write_log("Attempt to instantiate too many instances", GM_FATAL)
    end if
    if (associated(registered_models(num_models)%p)) then
       ! This most likely indicates a programming error, not a user error
       call write_log("Somehow we are attempting to register an already-registered model", &
            GM_FATAL)
    end if

    registered_models(num_models)%p => model
    model%model_id = num_models
  end subroutine

  subroutine check_num_models()
    !> Checks that the number of calls to register_model matched the expected number as
    !> set in the call to set_max_models. It only makes sense to call this if
    !> set_max_models was called (and assumes that the max set there is the expected
    !> actual number of models in this run); if this is called, it should be called after
    !> all calls to register_model.
    character(len=256) :: message

    if (num_models /= max_models) then
       write(message,'(a,i0,a,i0,a)') 'Number of calls to register_model (', num_models, &
            ') does not match max_models (', max_models, ')'
       call write_log(trim(message), GM_FATAL)
    end if
  end subroutine check_num_models

  subroutine deregister_model(model)
    !> Removes a model from the registry. Normally this should only be done when
    !> glide_finalise is called on the model, and is done automatically by that function.
    !>
    !> This does NOT decrement num_models, so weird behavior could arise if trying to use
    !> num_models, call register_model, etc. after calling deregister_model. As long as
    !> deregister_model is only called at the end of the run, this shouldn't be a problem.
    type(glide_global_type) :: model

    if (model%model_id < 1 .or. model%model_id > num_models) then
        call write_log("Attempting to deregister a non-allocated model", GM_WARNING)
    else
        registered_models(model%model_id)%p => null()
        model%model_id = 0
    end if
  end subroutine

  integer function get_num_models()
    !> Get the number of models actually registered so far.
    !>
    !> Typically this should only be called after all calls to register_model have been
    !> made.
    get_num_models = num_models
  end function get_num_models

  integer function get_max_models()
    !> Get the maximum number of models allowed in this run.
    !>
    !> If set_max_models is called at some point in initialization, get_max_models should
    !> be called after set_max_models (otherwise get_max_models will return the default,
    !> hard-coded max).
    get_max_models = max_models
  end function get_max_models

end module glide_model_registry
