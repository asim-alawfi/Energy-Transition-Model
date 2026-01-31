function cond=permfix(kind,name,dim,comp,varargin)
default={{'perm','trafo'},[],{'kron','node_dim'},1,'t',[0,0.1],'file',''};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
prm=options.perm(:);
nt=length(prm);
mat_int=@(ir,ic,nr,nc)kron(full(sparse(ir(:),ic(:),ones(numel(ir),1),nr,nc)),options.kron);
stateproj=mat_int(1:nt,prm,nt,dim);
proj=1:nt;
if isempty(prm)
    cond=feval(['dde_',kind,'_lincond'],name,dim*options.kron,comp,'trafo',zeros(0),'stateproj',stateproj);
    return
end
t_applications=struct('stst',{{}},'psol',{{'condprojint',options.t(:)*[1,1]}});
t_apply=t_applications.(kind);
trafomat=mat_int(proj([2:end,1]),proj,nt,nt);
cond=feval(['dde_',kind,'_lincond',options.file],name,dim*options.kron,comp,...
    'trafo',trafomat,'stateproj',stateproj,t_apply{:},pass_on{:});
end
