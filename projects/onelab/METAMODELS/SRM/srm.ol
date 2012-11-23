# Metamodel description
# 
Mesh.register(native, gmsh);
Mesh.in(OL.get(Arguments/FileName).geo);
Mesh.out(OL.get(Arguments/FileName).msh); 
Mesh.run(OL.get(Arguments/FileName).geo -v 1);
Mesh.frontPage(OL.get(Arguments/FileName).geo, OL.get(Arguments/FileName).msh);

Getdp/1ResolutionChoices.string(MagSta);
Getdp/1ResolutionChoices.setVisible(0);
Getdp/2PostOperationChoices.string(Torque);
Getdp/2PostOperationChoices.setVisible(0);

Getdp.register(native, getdp);
Getdp.in(OL.get(Arguments/FileName).pro );
Getdp.out(a.pos);
Getdp.run(OL.get(Arguments/FileName).pro -v 1);
Getdp.merge(a.pos);

