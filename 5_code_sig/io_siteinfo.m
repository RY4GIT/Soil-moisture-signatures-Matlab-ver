function [depth, lutype, nstation, npstation] = io_siteinfo(site)
    % A function that returns desired station information
    
    switch site
        case "WB" % site name abbreviation
            depth = [5; 20; 50]; % soil moisture sensor depth  
            nstation = 150; % number of total soil moisture stations 
            npstation = 1; % number of precipitation stations 
            lutype =  ["Forested";"Deforested"];  % land-use type

        case "HB"
            depth = [5; 10; 40; 80; 160]; 
            nstation = 10; 
            npstation = 1;
            lutype = ["Greenspace";"Urban"];

        case "RM"
            depth = [5; 10; 20; 40; 80]; 
            nstation = 15; 
            npstation = 1;
            lutype = ["Deep GW";"Shallow GW"];

        case "TX"
            depth = [5; 10; 20; 50]; 
            nstation = 40; 
            npstation = 40;
            lutype = ["Ungrazed";"Grazed"];

        case "MQ"
            depth = [5; 10; 20; 40; 80]; 
            nstation = 20; 
            npstation = 1;
            lutype = ["Non-wetland"; "Wetland"];


        case "OZ"
            depth = [3; 4; 15; 45; 75]; 
            nstation = 38;  
            npstation = 38;
            lutype = ["Grass";"Grazed";"Crop"];

    end

end