import stratiMesh

mesh = stratiMesh.stratiMesh(folder='../', xdmfName = 'stratal_series', dispTime=500.)
mesh.outputSteps(startTime=0.,endTime=50000.)

