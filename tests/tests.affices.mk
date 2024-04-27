affix.basic         := basic
affix.utility       := utility basic external
affix.combinatorics := combinatorics basic external
affix.operators     := operators combinatorics utility basic external
affix.dense         := dense sparse operators combinatorics utility basic external
affix.sparse        := sparse dense operators combinatorics utility basic external
affix.solver        := solver sparse dense operators combinatorics utility basic external
affix.mesh          := mesh dense sparse operators combinatorics utility basic external
affix.vtk           := vtk $(affix.mesh)
affix.fem           := fem mesh solver dense sparse operators combinatorics utility basic external
affix.solverfem     := fem mesh solver dense sparse operators combinatorics utility basic external
affix.sullivan2D    := fem vtk mesh solver dense sparse operators combinatorics utility basic external
affix.sullivan3D    := fem vtk mesh solver dense sparse operators combinatorics utility basic external
affix.whitney2D     := fem vtk mesh solver dense sparse operators combinatorics utility basic external
affix.whitney3D     := fem vtk mesh solver dense sparse operators combinatorics utility basic external