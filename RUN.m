% random walk on unweigted network
% the method uses the standard network diffusion algorithm based on Personalized PageRank (PPR): 
% Lawrence Page, Sergey Brin, Rajeev Motwani, and Terry Winograd. The pager-ank citation ranking: bringing order to the web. 1999
% Hyunghoon Cho, Bonnie Berger, and Jian Peng. Diffusion Component Analysis: Unraveling Functional Topology in Biological Networks. 2015
% The subrountine of implementing random walk was from: https://github.com/hhcho/diffusion-component-analysis/blob/master/code/run_diffusion.m

% A. It first computes reachability of a given node on the network to any other nodes on the network.
% B. It then labels CHD-associated proteins on the network. 
% C. For a given protein on the network, we asked if this protein is more reachable to these CHD proteins relative to all proteins on the network.
% D. We considered %significant genes/proteins with a fold change of reachability greater than 2 and corrected P value (Storey's FDR<=0.05).

% [Input]
% A: adjacency matrix (could be weighted)
% geneName: gene names of adjacency matrix
% chdGenes: known CHD genes list
% method: 'personalized-pagerank' or 'self-diffusion'
% maxiter: max number of iterations
% rsp: restart probability  
%
% [Output]
% SS0: the identified proteins on the network, excluding known CHD proteins

A = matrix;
method = 'personalized-pagerank';
geneName = geneName;
chdGenes = chdGenes;
maxiter = 500;
rsp = 0.1;

SS0 = RUN(A, method, geneName, chdGenes, maxiter, rsp);

function [SS] = RUN(A, method, geneName, chdGenes, maxiter, rsp)

% QA is the probablity matrix quantifying reachability of a given node to all other nodes on the network
QA = run_diffusion(A, method, struct('maxiter', maxiter, 'reset_prob', rsp));

% chdIdx is the indices of known CHD associated proteins on the network
[~,chdIdx]=intersect(geneName,chdGenes);

% Qchd is the reachability matrix for any node to a set of pre-defined known CHD associated proteins
Qchd=QA(:,chdIdx);

Q=(mean(QA,2)); % the average reachability of a given protein to all other proteins
Q0=(mean(Qchd,2)); % the average reachability of a given protein to CHD proteins

FC=Q0./Q; % fold chanage of reachability

%
P=[];
for i = 1:numel(geneName) % # of proteins on the network

    P(i) = ranksum(QA(i,:),Qchd(i,:)); % to test if a protein is more reachable to CHD proteins than background

end

fdr = mafdr(P)';

SS = setdiff(geneName(fdr<=0.05&FC>=2),chdGenes); %SS is the identified proteins on the network, excluding known CHD proteins

end
	
    
function [Q] = run_diffusion(A, method, aux)

    n = size(A, 1);
    renorm = @(M) bsxfun(@rdivide, M, sum(M));

    A = A + diag(sum(A) == 0); % Add self-edges to isolated nodes
    P = renorm(A);
    fprintf('rsp=%f\n',aux.reset_prob);
    if strcmp(method, 'personalized-pagerank')
      assert(isfield(aux, 'maxiter'))
      assert(isfield(aux, 'reset_prob'))
      
      reset = eye(n);
      Q = reset;
      for i = 1:aux.maxiter
        Q_new = aux.reset_prob * reset + (1 - aux.reset_prob) * P * Q;
        delta = norm(Q - Q_new, 'fro');
         fprintf('Iter %d. Frobenius norm: %f\n', i, delta);
        Q = Q_new;
        if delta < 1e-6
           fprintf('Converged.\n');
          break
        end
      end

    elseif strcmp(method, 'self-diffusion')
      assert(isfield(aux, 'maxiter'))

      Q = eye(n);
      for i = 1:aux.maxiter
        Q_new = eye(n) + P * Q;
         fprintf('Iter %d. Frobenius norm: %f\n', i, norm(renorm(Q) - renorm(Q_new), 'fro'));
        Q = Q_new;
      end

    else
      error(['Unknown method: ', method])
    end
    Q = bsxfun(@rdivide, Q, sum(Q));

end
