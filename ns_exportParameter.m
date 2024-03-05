%% Export parameter values to text files

function ns_exportParameter(gmp,mcp,odir,rpc)

fileId=fopen(sprintf('%sr%d_parameter.txt',odir,rpc),'w');

fprintf(fileId,sprintf('Box size(bs):%s\n',num2str(gmp.bs)));
fprintf(fileId,sprintf('Edge segment criterion(edpc):%s\n',num2str(gmp.edpc)));
fprintf(fileId,sprintf('Short edge criterion(shEd):%s\n\n',num2str(gmp.shEd)));

fprintf(fileId,sprintf('Run time(runTm in a unit of taut):%s\n',num2str(gmp.runTm)));
fprintf(fileId,sprintf('Simulation time step(dt):%s\n',num2str(gmp.dt)));
fprintf(fileId,sprintf('Data point time step(dt*svPt):%s\n\n',num2str(gmp.dt*gmp.svPt)));

fprintf(fileId,sprintf('Density parameter(psi):%s\n',num2str(mcp.psi)));
fprintf(fileId,sprintf('Nuclear packing parameter(nvf):%s\n\n',num2str(mcp.nvf)));
fprintf(fileId,sprintf('Nuclear aspect ratio(asp):%s\n\n',num2str(gmp.asp)));

fprintf(fileId,sprintf('Tension fluctuation times scale(taut):%s\n',num2str(mcp.taut)));
fprintf(fileId,sprintf('Nuclear translational relaxation times scale(taun):%s\n\n',num2str(mcp.taun)));
fprintf(fileId,sprintf('Nuclear rotational relaxation times scale(taua):%s\n\n',num2str(mcp.taua)));


fprintf(fileId,sprintf('Tension fluctuation(mu):%s\n',num2str(mcp.mu)));
fprintf(fileId,sprintf('Nuclear fluctuation(mun):%s\n\n',num2str(mcp.mun)));

fprintf(fileId,sprintf('Nuclear stiffness(nstf):%s\n',num2str(mcp.nstf)));
fprintf(fileId,sprintf('Nuclear restoring force(nrst):%s\n',num2str(mcp.nrst)));
fprintf(fileId,sprintf('Tension to pressure ratio(del):%s\n',num2str(mcp.del)));

fclose(fileId);

end