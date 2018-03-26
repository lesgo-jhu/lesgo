# Branches
The master branch is reserved for stable snapshots of the code that have been proven and compile and run for a large subset of configurations. Before committing changes to master, run the bundled test-lesgo script. The modified code must pass all compilation and runtime testing with no errors. Warnings at compile time are discouraged, but allowed.

Active development occurs in the "devel_branch" branch. For major code changes, consider creating a new branch first.

# Style guide
## License Heading
   All source files must begin with the licensing information:
    
    !!
    !!  Copyright (C) xxxx-xxxx  Johns Hopkins University
    !!
    !!  This file is part of lesgo.
    !!
    !!  lesgo is free software: you can redistribute it and/or modify
    !!  it under the terms of the GNU General Public License as published by
    !!  the Free Software Foundation, either version 3 of the License, or
    !!  (at your option) any later version.
    !!
    !!  lesgo is distributed in the hope that it will be useful,
    !!  but WITHOUT ANY WARRANTY; without even the implied warranty of
    !!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !!  GNU General Public License for more details.
    !!
    !!  You should have received a copy of the GNU General Public License
    !!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
    !!

## Maximum line width
All lines must be less than 80 characters in width. A long line in fortran
source files should be continued with an ampersand at column 80. For example:

    k = (1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10)                                   &
        * (1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10)
        
## Tabs vs. spaces
Use spaces for all alignments. Indentation uses 4 spaces instead of a tab.
Set your text editor to emit 4 spaces when you press the tab key. Do loops
should be indented 4 spaces:

    do i = 1, N
        {Do something}
    end do

Sequential loops may be aligned together:

    do i = 1, Nx
    do j = 1, Ny
        {Do something}
    end do
    end do

However, use this style if the do loops do not appear together:

    do i = 1, Nx
        {Do something}
        do j = 1, Ny
            {Do something else}
        end do
    end do

Do not indent inside module, subroutine, and function blocks

## Comments for modules and procedures.
Demarcate modules and precures using a ***'s and place the description of the
block directly below. The ***s should extend to the end of the column
(80 character width). For example:

    !*******************************************************************************
    function foo() result(bar)
    !*******************************************************************************
    !
    ! This function returns bar when calling foo. bar = foo()
    !

## Preprocessor directives
Use C preprocessor directives to turn features on and off at compile time.
In the Fortran implementation, the directives must not be indented. As a
result, indent other portions of the code as if these directives were not
there. For example:

    if (dothis = .true.) then
        write(*,*) "dothis is true"
    #ifdef PPFLAG
    else
        write(*,*) "HELLO FLAG!"
    #endif
    endif

If flag is defined, the preprocessor would write a file like this:

    if (dothis = .true.) then
        write(*,*) "dothis is true"
    else
        write(*,*) "HELLO FLAG!"
    endif

if PPFLAG is not defined, the preprocessor would write a file like this:

    if (dothis = .true.) then
        write(*,*) "dothis is true"
    endif

## Modules
All modules that need to have memory allocated or preliminary values
calculated must have a subroutine called MODULENAME_init that is called in
initialize().

Deallocation of memory must be performed in a subrotuine MODULENAME_finalize
that is called in finalize().
