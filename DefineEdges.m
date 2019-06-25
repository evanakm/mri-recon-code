% MRI BRAIN PHANTOM
% Generated with ExportPhantom.m on 02-May-2011 17:44:31
% Contact: Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne
Edges.FOV = [1,1];
Edges.region = cell(2*109,1);
%Brain.region = cell(2,1);
%epsilon = .03;
% Gray level 0: Black background
for i=1:109
    if i >= 11
        epsilon = .03;
    else
        epsilon = .03;
    end
    
    if(strcmp(Brain.region{i}.type,'bezier'))
        Edges.region{(2*i)-1} = struct('type','bezier','weight',1,'control',(1+epsilon)*Brain.region{i}.control);
        Edges.region{2*i} = struct('type','bezier','weight',-1,'control',(1-epsilon)*Brain.region{i}.control);
    else
        Edges.region{(2*i)-1} = struct('type','bezier','weight',1,'control',(1+.03)*Brain.region{1}.control);
        Edges.region{2*i} = struct('type','bezier','weight',-1,'control',(1-.03)*Brain.region{1}.control);
    end        
end

clear epsilon;