% This function will return the Hamiltonian for the Kitaev model on a
% surface of genus G+1 with a vortex configuration described by V. The
% surface is created by joining G octagonal pieces of lattice together.
% The size of each octagonal piece is determined by 3 numbers, Na - the
% number of sites along a vertical edge, Nb - the number of sites along a
% diagonal edge, Nc - the number of sites along a horizontal edge.


% To ensure correct checker board pattern, should have Na - even, 
% Nb - even/odd, Nc - odd. These numbers, along with the total number of
% sites per octagonal piece Ntot and the number of ocetagnal piece G are
% stored in N

% The homology sector of the system is determined by the eigenvalues of the
% loop operators which are stored in L as shown below.

%      La = L(end-1)            Lb = L(end)              Lc = L(2*g)              Ld = L(2*g-1)  
%      __|________              ___________              ___________              ___________  
%     // |       \\            //         \\            //         \\            //         \\
%    //  |        \\          //           \\          //|          \\          //          |\\
%   //   |         \\        //             \\        // |           \\        //           | \\
%  //    |          \\      //               \\      //  |            \\      //            |  \\
% ||     |           ||    ||_________________||    ||   |_____________||    ||_____________|   ||
% ||     |           ||    ||                 ||    ||                |||    |||                ||
% ||     |           ||    ||                 ||    ||                |||    |||                ||
% ||     |           ||    ||                 ||    ||                |||    |||                ||
%  \\    |          //      \\               //      \\               //      \\               //
%   \\   |         //        \\             //        \\             //        \\             //
%    \\  |        //          \\           //          \\           //          \\           //
%     \\_|_______//            \\_________//            \\_________//            \\_________// 



function H =  Hamiltonian(Jx,Jy,Jz,L,k,V,N)

    Na = N(1); Nb = N(2); Nc = N(3); Ntot = N(4); G = N(5);
    
    delta_x = spalloc(G*Ntot,G*Ntot,G*Ntot);
    delta_y = spalloc(G*Ntot,G*Ntot,G*Ntot);
    delta_z = speye(G*Ntot,G*Ntot);
    
    % turn on the x and y links
    for g = 1:G
        
        for x=0:(2*Nb+Nc)-1
            if ( x<Nb-1 ) || ( x>=Nb+Nc && x< 2*Nb+Nc-1 )
                for y=0:Nb+Na-1
                    delta_x(IndxGTorus(x,y,g,N),IndxGTorus(x+1,y,g,N))=1;
                    delta_y(IndxGTorus(x,y,g,N),IndxGTorus(x,y+1,g,N))=1;
                end
            elseif x>=Nb && x<Nb+Nc-1
                for y=0:2*Nb+Na-1
                    delta_x(IndxGTorus(x,y,g,N),IndxGTorus(x+1,y,g,N))=1;
                    delta_y(IndxGTorus(x,y,g,N),IndxGTorus(x,y+1,g,N))=1;
                end
            elseif x == Nb-1
                for y=0:Nb+Na-1
                    delta_y(IndxGTorus(x,y,g,N),IndxGTorus(x,y+1,g,N))=1;
                    if y<=Nb
                        delta_x(IndxGTorus(x,y,g,N),IndxGTorus(x+1,y,g,N))=1;
                    else
                        delta_x(IndxGTorus(x,y,g,N),IndxGTorus(Nb+Nc,y,g+1,N))=1;
                    end
                end
            elseif x == Nb+Nc-1
                for y=0:2*Nb+Na-1
                    delta_y(IndxGTorus(x,y,g,N),IndxGTorus(x,y+1,g,N))=1;
                    if y<=Nb
                        delta_x(IndxGTorus(x,y,g,N),IndxGTorus(x+1,y,g,N))=1;
                    else
                        delta_x(IndxGTorus(x,y,g,N),IndxGTorus(0,y-Nb,g,N))=1;
                    end
                end
            elseif x == 2*Nb+Nc-1
                for y=0:Nb+Na-1
                    delta_y(IndxGTorus(x,y,g,N),IndxGTorus(x,y+1,g,N))=1;
                    if y<1
                        delta_x(IndxGTorus(x,y,g,N),IndxGTorus(x+1,y,g,N))=1;
                    else
                        delta_x(IndxGTorus(x,y,g,N),IndxGTorus(Nb,Nb+y,g,N))=1;
                    end
                end
            end
        end
        
    end
    
    
    
    M = size(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Closed loops
    for g = 1:G
        

%       first corner
%       __|_________                ____________               __|_________        
%      //:|         \\             //::::::::::\\             // |::::::::\\
%     //|:|         |:\           //::::::::::::\\           //| |:::::::::\\
%    // |:|         |::\         //::::::::::::::\\         //:| |::::::::::\\
%   //  |:|         |::|\       //::::::::::::::::\\       //::| |:::::::::::\\
%  //   |:|         |::|\\     //::::::::::::::::::\\     //:::| |:::::::::::|\\
% ||  __|:|         |::| ||   ||::::::::::::::::::::||   ||::::|_|:::::::::::| ||
% || |::::|              ||   ||::::::::::::::::::::||   ||::::::::::::::::::| ||
% || |::::|              ||   ||::::::::::::::::::::||   ||::::::::::::::::::| ||
% || |::::|              ||   ||::::::::::::::::::::||   ||::::::::::::::::::| ||
%  \\|::::|             //     \\:::::::::::::::::://     \\:::|-|:::::::::::|//
%   \\::::|            //       \\:::::::::::::::://       \\::| |::::::::::://
%    \\:::|           //         \\:::::::::::::://         \\:| |:::::::::://
%     \\::|          //           \\:::::::::::://           \\| |::::::::://
%      \\:|_________//             \\:::::::::://             \\_|:::::::://        
        
        % first corner
        %%
        if g>1
            
            v = 0;
            for h=1:M(1)
                    
                if (V(h,1)<Nb-1 || (V(h,1)>=Nb+Nc && V(h,1)<2*Nb+Nc-1) || (V(h,1)>=Nb && V(h,1)<=Nb+Nc-1)) && V(h,3)==g
                    v = v+1;
                        
                elseif ((V(h,1)>=Nb+Nc && V(h,1)<2*Nb+Nc-1) || V(h,1)==Nb-1 || V(h,1)==2*Nb+Nc-1) && V(h,3)==1
                    v = v+1;
                        
                elseif V(h,3)<g && V(h,3)>1
                    v = v+1;
                end
                    
            end
            % correct for over counting vortices
            for j= Nb+Na-1:-1:1
                for h=1:M(1)
                    
                    if (V(h,1)==2*Nb+Nc-1 && V(h,2)==j) && V(h,3) == g
                        v = v+1;
                    end
                    
                end
                delta_x(IndxGTorus(2*Nb+Nc-1,j,g,N),:) = (-1)^(Nb+Nc)*((-1)^v)*L(end-1)*L(2*g)*L(1)*L(2)*delta_x(IndxGTorus(2*Nb+Nc-1,j,g,N),:);
            end
        else
            v = 0;
            for h=1:M(1)
                if (V(h,1)>=Nb+Nc && V(h,1)< 2*Nb+Nc-1) && V(h,3) == 1
                    v = v+1;
                end
            end
                
            for j= Nb+Na-1:-1:1
                for h=1:M(1)
                    if (V(h,1)==2*Nb+Nc-1 && V(h,2)==j) && V(h,3) == 1
                        v = v+1;
                    end
                end
                delta_x(IndxGTorus(2*Nb+Nc-1,j,1,N),:) = (-1)^(Nb+Nc)*((-1)^v)*L(end-1)*L(1)*delta_x(IndxGTorus(2*Nb+Nc-1,j,1,N),:);
            end
        end
        %%
        
%       second corner
%       __|_________                ____________               ________|___        
%      //:|         \\             //::::::::::\\             //       |::\\
%     //|:|         |:\           //::::::::::::\\           //        |::|\\
%    // |:|         |::\         //::::::::::::::\\         //         |::| \\
%   //  |:|         |::|\       //::::::::::::::::\\       //          |::|  \\
%  //   |:|         |::|\\     //::::::::::::::::::\\     //           |::|   \\
% ||  __|:|         |::| ||   ||::::::::::::::::::::||   ||            |::|    ||
% || |::::|              ||   ||::::::::::::::::::::||   ||                    ||
% || |::::|              ||   ||::::::::::::::::::::||   ||                    ||
% || |::::|              ||   ||::::::::::::::::::::||   ||             __     ||
%  \\|::::|             //     \\:::::::::::::::::://     \\           |::|   //
%   \\::::|            //       \\:::::::::::::::://       \\          |::|  //
%    \\:::|           //         \\:::::::::::::://         \\         |::| //
%     \\::|          //           \\:::::::::::://           \\        |::|//
%      \\:|_________//             \\:::::::::://             \\_______|::// 
        
        % second corner
        %%
        if g>1
            
            v = 0;
            for h=1:M(1)
                    
                if (V(h,1)==Nb+Nc-1 && V(h,2)<=Nb) && V(h,3)==g
                    v = v+1;
                        
                elseif ((V(h,1)>=Nb+Nc && V(h,1)<2*Nb+Nc-1) || V(h,1)==Nb-1 || V(h,1)==2*Nb+Nc-1 ) && V(h,3)==1
                    v = v+1;
                        
                elseif V(h,3)<g && V(h,3)>1
                    v = v+1;
                end
                    
            end
            
            for j= Nb:2*Nb+Na-2
                for h=1:M(1)
                    
                    if (V(h,1)==Nb+Nc-1 && V(h,2)==j) && V(h,3) == g
                        v = v+1;
                    end
                    
                end
                delta_x(IndxGTorus(Nb+Nc-1,j+1,g,N),:) = (-1)^(Nb+Nc)*((-1)^v)*L(end-1)*L(2*g-1)*L(1)*L(2)*delta_x(IndxGTorus(Nb+Nc-1,j+1,g,N),:);
            end
            
        else
            v = 0;
            for h=1:M(1)
                if (V(h,1)<Nb+Nc-1 && V(h,1) ~= Nb-1) && V(h,3) == 1
                    v = v+1;
                end
            end   
            for j= 2*Nb+Na-1:-1:Nb+1
                for h=1:M(1)
                    if (V(h,1)==Nb+Nc-1 && V(h,2)==j) && V(h,3) == 1
                        v = v+1;
                    end
                end
                delta_x(IndxGTorus(Nb+Nc-1,j,1,N),:) = (-1)^(Nb+Nc)*((-1)^v)*L(end-1)*L(2)*delta_x(IndxGTorus(Nb+Nc-1,j,1,N),:);
            end
        end
        %%
        
%       third corner
%       ____________                ____________               ____________        
%      //           \\             //          \\             //          \\
%     //             \\           //            \\           //            \\
%    //               \\         //              \\         //              \\
%   //                 \\       //                \\       //                \\
%  //                   \\     //                  \\     //                  \\
% ||                     ||   ||                   _||   ||_                   ||
% ||                     ||   ||                  |:||   ||:|                  ||
% ||                     ||   ||                  |:||   ||:|                  ||
% ||                     ||   ||                  |_||   ||_|                  ||
%  \\                   //     \\                  //     \\                  //
%   \\                 //       \\                //       \\                //
%    \\               //         \\              //         \\              //
%     \\             //           \\            //           \\            //
%      \\___________//             \\__________//             \\__________// 
        
        % third corner
        %%
        if g < G
            for j= Nb+Na-1:-1:Nb+1
                for h=1:M(1)
                    if (V(h,1)==Nb-1 && V(h,2)==j) && V(h,3) == 1
                        v = v+1;
                    end
                end
                delta_x(IndxGTorus(Nb-1,j,g,N),:) = ((-1)^v)*L(2*g)*L(2*(g+1)-1)*delta_x(IndxGTorus(Nb-1,j,g,N),:);
            end
            
        else
            
            for j= Nb+Na-1:-1:Nb+1
                for h=1:M(1)
                    if (V(h,1)==Nb-1 && V(h,2)==j) && V(h,3) == G
                        v = v+1; 
                    end
                end
                delta_x(IndxGTorus(Nb-1,j,G,N),:) = (-1)^(Nc)*((-1)^v)*L(1)*L(2*G)*L(end)*delta_x(IndxGTorus(Nb-1,j,G,N),:);
            end
            delta_x(IndxGTorus(2*Nb+Nc-1,0,G,N),:) = (-1)^(Nc)*L(end)*delta_x(IndxGTorus(2*Nb+Nc-1,0,G,N),:);
            
        end
        

        %%
        
        
%           (B)                         (A)                        (C)        
%       __|__________               ____________               ____________        
%      // |::::|    \\             //          \\             //          \\
%     //  |::::|     \\           /:|           \\           //           |:\
%    //   |::::|      \\         /::|            \\         //            |::\
%   //    |::::|       \\       /|::|             \\       //             |::|\
%  //     |::::|        \\     //|::|              \\     //              |::|\\
% ||      |::::|         ||   || |::|______________ ||   || ______________|::| ||
% ||      |::::|         ||   ||               |::| ||   || |::|               ||
% ||      |::::|         ||   ||               |::| ||   || |::|               ||
% ||      |::::|         ||   ||               |::| ||   || |::|               ||
%  \\     |::::|        //     \\              |::|//     \\|::|              //
%   \\    |::::|       //       \\             |:://       \\::|             //
%    \\   |::::|      //         \\            |://         \\:|            //
%     \\  |::::|     //           \\           |//           \\|           //
%      \\_|::::|____//             \\__________//             \\__________// 
        
               
        % y-links (A)
        %%
        for x=0:Nb-1
            v = 0;
            for h=1:M(1)
                if ((V(h,1)<Nb-1) && V(h,1)>=x) && V(h,3) == g
                    v = v+1;
                end
            end
            delta_y(IndxGTorus(x,Nb+Na-1,g,N),:) = (-1)^(Nb+Nc)*((-1)^v)*L(2*g)*delta_y(IndxGTorus(x,Nb+Na-1,g,N),:);
        end
        %%
    
    
        % y-links (B)
        %%
        if g>1
            for x=Nb:Nb+Nc-1
                v = 0;
                
                for h=1:M(1)
                    if ( (( (V(h,1)>=x) && V(h,1)<Nb+Nc-1) || V(h,1)<Nb-1 || V(h,1)==Nb+Nc-1 ) && V(h,3)==g ) || (V(h,3)<g && V(h,3)>1)
                        v = v+1;
                    elseif ((V(h,1)>=Nb+Nc && V(h,1)<2*Nb+Nc-1) || V(h,1)==Nb-1 || V(h,1)==2*Nb+Nc-1) && V(h,3)==1
                        v = v+1;
                    end
                end
                
                delta_y(IndxGTorus(x,2*Nb+Na-1,g,N),:) = ((-1)^v)*L(end-1)*L(2*g)*L(2*g-1)*L(1)*L(2)*delta_y(IndxGTorus(x,2*Nb+Na-1,g,N),:);
            end
        else
            for x=Nb:Nb+Nc-1
                v = 0;
                for h=1:M(1)
                    if ( ((V(h,1)<x) && V(h,1)>=Nb) && V(h,3) == 1 )
                        v = v+1;
                    end
                end
                delta_y(IndxGTorus(x,2*Nb+Na-1,g,N),:) = ((-1)^v)*L(end-1)*delta_y(IndxGTorus(x,2*Nb+Na-1,g,N),:);
            end
        end
        %%
    
    
        % y-links (C)
        %%
        for x=Nb+Nc:2*Nb+Nc-1
            v = 0;
            for h=1:M(1)
                if ((V(h,1)>=Nb+Nc) && V(h,1)<x) && V(h,3) == g
                    v = v+1;
                end
            end
            delta_y(IndxGTorus(x,Nb+Na-1,g,N),:) = (-1)^(Nb+Nc)*((-1)^v)*L(2*g-1)*delta_y(IndxGTorus(x,Nb+Na-1,g,N),:);
        end        
        
    end
    
    
    % flip the sign of the links appropriately to encode the the values of
    % the X operators, according to the vortex configuration V
    for h=1:M(1)
        
        if (V(h,1) < Nb-1) || ( V(h,1)>=Nb+Nc && V(h,1)<2*Nb+Nc-1 )
            for j=floor(V(h,2))+1:Nb+Na-1
                delta_x(IndxGTorus(floor(V(h,1)),j,V(h,3),N),:) = -1*delta_x(IndxGTorus(floor(V(h,1)),j,V(h,3),N),:);
            end
        elseif ( V(h,1) >= Nb ) && ( V(h,1)<Nb+Nc-1 )
            for j=floor(V(h,2))+1:2*Nb+Na-1
                delta_x(IndxGTorus(floor(V(h,1)),j,V(h,3),N),:) = -1*delta_x(IndxGTorus(floor(V(h,1)),j,V(h,3),N),:);
            end
        elseif (V(h,1)==Nb-1 || V(h,1)==Nb+Nc-1) && (V(h,2)<Nb)
            for j=floor(V(h,2))+1:Nb
                delta_x(IndxGTorus(floor(V(h,1)),j,V(h,3),N),:) = -1*delta_x(IndxGTorus(floor(V(h,1)),j,V(h,3),N),:);
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    Xi = 2*Jz*delta_z + Jx*( delta_x'+delta_x ) + Jy*(delta_y'+delta_y);
    Del = Jx*(delta_x'-delta_x) + Jy*(delta_y'-delta_y);
    
    % create the three body magnetic terms
    %     P(1)                     P(2)                    P(3)                    P(4)
    del = -1i*(delta_x-delta_x') + 1i*(delta_y-delta_y') - 1i*(delta_x-delta_x') + 1i*(delta_y-delta_y');
    
    %         P(5)                                            P(6)
    xi  =     - 1i*(delta_x*delta_y'- (delta_x*delta_y')' ) + 1i*( delta_y'*delta_x - (delta_y'*delta_x)' );
    del = del + 1i*(delta_x*delta_y'- (delta_x*delta_y')' ) + 1i*( delta_y'*delta_x - (delta_y'*delta_x)' );
    
    % include the magnetic terms in the BdG blocks
    Xi = Xi + k.*xi;
    Del = Del + k.*del;
    
    H = 0.5*[ Xi,Del ; Del',-(Xi).' ];

end