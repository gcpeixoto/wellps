function dv = setNeighDist(strel)
%SETNEIGHDIST chooses the structuring element that
%             drives the neighbourhood connectivity.
%
% PARAMETERS: 
%       - strel: a char indicating the structuring element.
%                It accepts two values: '6' and '26'
%
% RETURNS: 
%       - dv: distance to use as criterion for connectivity
%
% 
% REMARK: '6' neighbours is the only rule truly acceptable to 
%         give dynamic connectivity in numerical simulators. 
%         At least, we believe that, just as for IMEX, other 
%         softwares also are based on finite volume or similar 
%         methods that take into account face connectivy to 
%         compute transmissibilities, fluxes, etc. 
%         The option for '26' neighbours is left here for 
%         mere utility. 
%
%{
     Developed at LaMEP/UFPB, Brazil
     @gcpeixoto

%}

% checking
assert(ischar(strel),'strel must be a char: %s or %s neighbours.','6','26');

switch strel
    
    % 6 neighbours (3D von Neumman's neighborhood) => distance = 1
    %       _
    %     _|_|_    |
    %    |_|_|_| - o -
    %      |_|     |
    %       
    %   i-1 i i+1
    %
    case '6'       
        dv = 1;        
    
        
    % 26 neighbours (3D Moore's neighborhood) => distance = sqrt(3)
    %     _ _ _
    %    |_|_|_|  \ | /
    %    |_|_|_|  - o -
    %    |_|_|_|  / | \
    %       
    %   i-1 i i+1
    %
    case '26'
        dv = sqrt(3);
        
end
        

end

