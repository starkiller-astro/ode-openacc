# Overview
`ode-openacc` is intended to be a self-contained set of Fortran source code
implementing basic first-order ordinary differential equation (ODE) integration
accelerated with OpenACC.  This allows for easy experimentation and testing of
ideas, which can inform development of production codes like
[Maestro](https://github.com/BoxLib-Codes/MAESTRO).  It can also be useful as a
so-called mini-app to be used for the verification of compilers implementing the
OpenACC standard.

# Requirements
The code is developed mostly on Fedora Linux machines and using PGI compilers,
so should work in similar environments.  If you have problems with the code,
please file an issue describing the issue.  We will also consider pull requests
fixing issues or adding functionality.

# A basic example
First, clone the repo.  At a command prompt, execute
```bash
git clone https://github.com/adam-m-jcbs/ode-openacc.git
```

Now build an executable.  The simplest example is in `first-order`.
```bash
cd ode-openacc/first-order
make
```

By default, this will build a version of the code that runs in serial on the
CPU using PGI compilers, giving results like
```bash
./test_react.pgf95_local.exe 
            0
   6.4935629937149325E-002   0.5000000000000000        0.4350643700628506     
 runtime:     11.69168210029602       seconds
```

To build with OpenACC, do
```bash
make ACC=t
./test_react.pgf95_local.acc.exe
```

Hopefully by examining this example and the relevant source/make infrastructure,
you can get a feel for the code and how to modify it for your purposes.
