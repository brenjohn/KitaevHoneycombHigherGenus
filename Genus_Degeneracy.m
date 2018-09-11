clear; 

% Set lattice dimensions
Na = 6; Nb = 6; Nc = 6;
Ntot = 2*Nb*(Nb+Na)+Nc*(2*Nb+Na);
G = 1;

% V holds the coorendinates of the vortices in the system
% For the vortex free sector V is an empty array 
V = [];

% Concatenate the dimensions of the lattice into a single variable for the 
% Hamiltonian function
N = [ Na, Nb, Nc, Ntot, G];

% Set model parameters
Ji = 0.0; Jf = 1.0; Jstep = 0.05;
k = 0.2; Jz = 1; 
J_len = length(Ji:Jstep:Jf); L_len = 2^(2*G+2);

% allocate memory to save ground state energies and which homology sectors
% become excited
raised = zeros(J_len,L_len);
En = zeros(G*Ntot,J_len,L_len);
GS_E = zeros(J_len,L_len);


% loop through all the homology sectors for each value of Ji:Jstep:Jf
for i = 1:J_len*L_len

    display( strcat('iteration_',num2str(i),'_out_of_',num2str(J_len*L_len)) )
    
    j = ceil(i/L_len);
    l = j*L_len-i;
    
    % Calculate the value of Jx=Jy=J and the eigenvalues for the Loop
    % operators for this particular homology sector
    J = Ji+Jstep*(j-1);
    L = 2*mod(dec2bin(l,2*(G+1)),2)-1;
    
    % Calculate the Hamiltonian
    H = Hamiltonian(J,J,Jz,L,k,V,N);
    [UV,E] = eig(full(H));
    
    % rearrange to have positive eigenvalues in increaseing order
    % and then negative
    E = diag( E([G*Ntot+1:2*G*Ntot,G*Ntot:-1:1],[G*Ntot+1:2*G*Ntot,G*Ntot:-1:1]) );
    UV = [UV(:,G*Ntot+1:2*N),UV(:,G*Ntot:-1:1)];
    
    % singular value decomposition of upper left block
    U = UV(1:G*Ntot,1:G*Ntot);
    [LU,DU,RU] = svd(U);
    
    du = diag(DU);
    II = find(abs(du)<0.000001);
    occupied = length(II);
    s = (mod(occupied,2)==1);
    
    % Calculate the parity dependence on loop operators
    % Note the absence of the (-1)^(Na*Nc) factor from the formula is due
    % to the fact that for odd Nc, L_1 [or L(end)] has an odd number of -1
    % factors. This implies the eigenvalues are flipped. ie for even Nc,
    % the eigenvalues of L_1 are +/- 1 but for odd Nc they are -/+1,
    % meaning the name of some homolgy sectors may change between even and
    % odd Nc
    ll = 1;
    for jj = 3:length(L)
        ll = ll*L(jj-2);
    end
    LoopParity = (-L(end))^Na*(-L(end-1)*L(1)*L(2))^(Nc*G)*(ll^Nc);
    
    % modify fermion parity according to homology sector
    if LoopParity==-1
        
        s = mod(s+1,2);
    end
    
    if s == 0
        GS_E( j,j*L_len-i+1 ) = -sum( E(1:G*Ntot) )/2;
    else
        GS_E( j,j*L_len-i+1 ) = -sum( E(1:G*Ntot) )/2 + E(1);
        raised(j,j*L_len-i+1) = 1;
    end
    En(:,j,j*L_len-i+1) = E(1:G*Ntot);
     
end

% Calculate the degeneracy
Degeneracy = length(raised(end,:))-length(find(raised(end,:)));

% Plot E-E_min vs J where E is the energy of the ground state in each
% homology sector
[m,I] = min(GS_E(floor(end/2),:));
plot(Ji:Jstep:Jf,GS_E-GS_E(:,I)*ones(1,4^(G+1)),'Color',[0    0.4470    0.7410],'LineWidth', 2)
hold


% Plot E-E_min vs J where E is the energy of the first excited state in each
% homology sector
E1 = zeros(J_len,L_len);
for j = 1:J_len
    for l = 1:L_len
        if raised(j,l)==0
            E1(j,l) = -sum( En(:,j,l) )/2 + En(1,j,l) + En(2,j,l);
        else
            E1(j,l) = -sum( En(:,j,l) )/2 + En(2,j,l);
        end
    end
end

plot(Ji:Jstep:Jf,E1-GS_E(:,I)*ones(1,4^(G+1)),'Color',[0    0.4470    0.7410],'LineWidth', 2)

ylim([-0.05 2])
set(gca, 'FontSize', 20)
xlabel('J','FontSize',20)
ylabel('E-E_{min}','FontSize',20)