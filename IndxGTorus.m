% This function caulculates the index of of a particular site of the
% lattice given its' x,y coordinates on the g-th copy of an octagonal piece
% of square lattice


function i = IndxGTorus( x,y,g,N )

Na = N(1); Nb = N(2); Nc = N(3); Ntot = N(4); G = N(5);

if g > G
    g = 1;
end

% G-torus boundary conditions
if x>(2*Nb+Nc)-1
    x = x-(2*Nb+Nc);
    g = mod(g,G)+1;
end

if x <= Nb-1
    
    if y>Nb+Na-1
        y = y-(Nb+Na);
        x = x+(Nb+Nc);
    end
    
    if x<=Nb-1
        i = y + (Nb+Na)*x;
    else
        % y should be 0 in this case
        i = (Nb+Na)*(x-(Nb+Nc)) + (Nb+Na)*Nb + (2*Nb+Na)*Nc;
    end
    
elseif x >= Nb && x <= Nb+Nc-1
    
    if y>(2*Nb+Na)-1
        y = y-(2*Nb+Na);
    end
    i = y + (2*Nb+Na)*(x-Nb) + Nb*(Nb+Na);
    
elseif x>=(Nb+Nc) && x<(2*Nb+Nc)
    
    if y>(Nb+Na-1)
        y = y-(Nb+Na);
        x = x-(Nb+Nc);%mod(x+(Nb+1),(2*Nb+Nc));
    end
    
    if x>=(Nb+Nc) && x<(2*Nb+Nc)
        i = y + (Nb+Na)*(x-(Nb+Nc)) + (Nb+Na)*Nb + (2*Nb+Na)*Nc;
    else
        i = (Nb+Na)*x;
    end
    
end

i = i+1+(g-1)*Ntot;

end