%%
clearvars
% close all
clc

% This script outputs data of Track distance, Track displacement, ...
% Track distance/displacement ratio, Step delta x & delta y, ...
% and Step velocity

%Examples see my thesis Fig3.7B,C,D

%Trackmate Plugin in Fiji, save the files as .csv to read them
%customise file names

%INPUT file name (one repeat at a time)
links = '190321_ctrl_links.csv'
spots = '190321_ctrl_spots.csv'

%OUTPUT results of this experiment (tip, vSIV, & branch)
Output_filename = '190321_ctrl_analysis.xlsx';
%Note there will be 3 empty defult sheets in the Excel file ...
% if using MATLAB verision <2019b

data=xlsread(links);
% read links file
data2=xlsread(spots);
% read spots file

% store information into arrays
track_ids = data(:,1);     %all cell tracks
cellID_source = data(:,2); %cell at ti
cellID_target = data(:,3); %cell at ti+1
colour = data(:,11);       %to sort tip, vSIV & branch cells
spotid = data2(:,1);       %individual cell ID
xcoord = data2(:,4);       % x coordinates 
ycoord = data2(:,5);       % y coordinates 
frames = data2(:,8);       % time point

%count number of tracks
uniq_track_ids = unique(track_ids);
length(track_ids)
numel(uniq_track_ids)
% 

dstn_m_all = [];         %tip cell  %distance of each track
dstn_b_all = [];         %vSIV cell
dstn_o_all = [];         %branch cell

dspl_m_all = [];                    %displacement of each track
dspl_b_all = [];
dspl_o_all = [];

ddr_m_all = [];                     %distance/displacement ratio
ddr_b_all = [];
ddr_o_all = [];

Tdelx_m_all = [];                     %Track deltax of tip cell
Tdely_m_all = []; 
Tdelx_b_all = [];                     %Track deltax of vSIV cell
Tdely_b_all = [];
Tdelx_o_all = [];                     %Track deltax of branch cell
Tdely_o_all = [];

delx_m =[];                   %delx of each step in the track
dely_m =[];
delx_b =[];
dely_b =[];                   %dely of each step in the track
delx_o =[];
dely_o =[];

sm =[];           %Velocity of tip cell
sb =[];
so =[];

%set some counters/flags to zero
proliferation_counter = 0; %count cell divisions
track_counter = 0; %count all sub_tracks
next_target_flag = 0;
i = 1;


for o = 1:numel(uniq_track_ids) %loop over each main track
    first_row = 1; %first row of main track is 'true'
    j=1; %reset subtrack number
    curr_trackid = track_ids(i);
    
    while track_ids(i) == curr_trackid %loop until we reach a new main track
        if first_row
            %this is the first row of this main track
             
            current_source = cellID_source(i);        %find spot id (from 'links.csv')
            k = find(spotid == cellID_source(i));     %in 'spot.csv'
            
            init_x_pos(j) = xcoord(k); %store initial x/y position for THIS subtrack
            init_y_pos(j) = ycoord(k);
            
                                  
            prev_x_pos(j) = xcoord(k); 
            prev_y_pos(j) = ycoord(k);
            
            t0 = frames(k);
            current_target = cellID_target(i);
            
            delx=[];
            dely=[];
            
            speed=[];
                        
            first_row = 0; %from now on this is false, until we reach a new track
            i=i+1;
            
        else
            if cellID_source(i) == current_target
                % if source (cell at ti+1) equals to target (cell at ti)
                % means it's the same cell at consecutive time points
                k = find(spotid == cellID_source(i));
                % then find the ID of the cell at ti+1
                
                tempdelx = xcoord(k) - prev_x_pos;
                tempdely = prev_y_pos - ycoord(k); %convert delty to ventral direction
                
                temp_dspl = sqrt(tempdelx.^2 + tempdely.^2);
                temp_speed = temp_dspl / 0.5 %30mins interval, turn unit into um/hour, speed of each step
                speed = [speed;temp_speed];
                
                delx = [delx;tempdelx]; %store delx
                dely = [dely;tempdely]; %store dely 
                
                              
                %then set new target
                prev_x_pos = xcoord(k); %(cell at ti+2)
                prev_y_pos = ycoord(k);
                
                current_target = cellID_target(i);
                current_i = i;
                i=i+1;
                
            else
                if i>1 && cellID_source(i) == cellID_source(i-1)
                    %a proliferation event has occurred at this position
                    % we want to mark this row as the start of a new
                    % subtrack to return to later
                    
                    proliferation_counter = proliferation_counter + 1;
                    
                    prev_dspl = sqrt(delx.^2 + dely.^2) %step displacement 
                    temp_dstn = sum(prev_dspl); %previous distance
                                                                                                  
                    if next_target_flag == 0    %no subtrack previously
                        
                        prev_dstn = temp_dstn; %store previous distance 
                    
                        prev_speed = speed;
                                                                       
                        next_target = cellID_target(i); %this will start the next sub_track to follow
                        next_target_row = i;
                        next_target_flag = 1;  %flag the first subtrack
                        next_source = cellID_source(i);
                        i=i+1;
                    else
                        %when a cell divides more than once, 
                        %but this rarely happens in the SIVP cells
                        %we already have a queue of events to return to!
                        
                        prev_dstn = [prev_dstn;temp_dstn]; %store previous distance 
                        
                        prev_speed = {prev_speed;speed}; %more than one array of speed
                                               
                        next_target = [next_target;cellID_target(i)]; %so append list
                        next_source = [next_source;cellID_source(i)];
                        next_target_row = [next_target_row;i]; %and append this list
                        next_target_flag = 1;
                        i=i+1;
                    end
                else
                    %ignore this row as it belongs to a different subtrack
                    i=i+1;
                end
            end
            
        end
        
        if i > length(track_ids)
            break
        end
        
    end
    
    %calculate trajectory properties of the first subtrack
    %or of the track of the cell if it does not divide
    
    k = find(spotid == current_target); %ID of the last cell
    tf = frames(k);
    tempdelx = xcoord(k) - prev_x_pos; 
    tempdely = prev_y_pos - ycoord(k); %convert delty to ventral direction
    
    delx = [delx;tempdelx]; %step delx
    dely = [dely;tempdely]; %step dely
    
    temp_dspl = sqrt(tempdelx.^2 + tempdely.^2);
    temp_speed = temp_dspl / 0.5 %30mins interval, turn unit into um/hour, speed of each step
    speed = [speed;temp_speed];
                 
    temp_track_dstn = sqrt(delx.^2 + dely.^2); %step displacement
    track_dstn = sum(temp_track_dstn); %track distance
                
    track_delx = xcoord(k) - init_x_pos %track delx
    track_dely = init_y_pos - ycoord(k)  %track dely
    
   
    track_dspl = sqrt(track_delx.^2 + track_dely.^2)  %track displacement
    dstn_dspl_ratio = track_dstn / track_dspl  %track distance/displacement ratio
   
%     sort results according to colour
    if colour(current_i)== -65281 %tip cell in magenta
                        
        ddr_m_all = [ddr_m_all;dstn_dspl_ratio];
        dstn_m_all = [dstn_m_all;track_dstn];
        dspl_m_all = [dspl_m_all;track_dspl];
        
        delx_m =[delx_m;delx];
        dely_m =[dely_m;dely];
        
        Tdelx_m_all = [Tdelx_m_all;track_delx];                     
        Tdely_m_all = [Tdely_m_all;track_dely]; 
        
        sm = [sm,speed]
        speed = [];
        
        delx=[];
        dely=[];
      
    else if colour(current_i) == -16776961 %vSIV cell in blue
        
        ddr_b_all = [ddr_b_all;dstn_dspl_ratio];
        dstn_b_all = [dstn_b_all;track_dstn];
        dspl_b_all = [dspl_b_all;track_dspl];

        delx_b =[delx_b;delx];
        dely_b =[dely_b;dely];
        
        Tdelx_b_all = [Tdelx_b_all;track_delx];                     
        Tdely_b_all = [Tdely_b_all;track_dely]; 
        
        sb = [sb,speed]
        speed = [];
        
        delx=[];
        dely=[];
        
    else colour(current_i) == -26368 %branch cell in orange
       
        ddr_o_all = [ddr_o_all;dstn_dspl_ratio];
        dstn_o_all = [dstn_o_all;track_dstn];
        dspl_o_all = [dspl_o_all;track_dspl];
        
        delx_o =[delx_o;delx];
        dely_o =[dely_o;dely];
        
        Tdelx_o_all = [Tdelx_o_all;track_delx];                     
        Tdely_o_all = [Tdely_o_all;track_dely]; 
        
        so = [so,speed]
        speed = [];
        
        delx=[];
        dely=[];
      
        end
    end
    
    %check to see if we had any proliferation events during that track, then return to search again from there.

    while next_target_flag
        if length(next_target_row) == 1
            %if we only have one proliferation event to deal with
            i=next_target_row;
            
            current_target = next_target; %start from the daughter cell
            current_source = next_source;
            
            ks = find(spotid == current_source);
            prev_x_pos = xcoord(ks);
            prev_y_pos = ycoord(ks);

             prev_speed_1 = prev_speed;
            ifc = iscell(prev_speed_1);
            if ifc == 1
                prev_speed_1 = cell2mat(prev_speed_1);
            end     
                       
            prev_dstn_1 = prev_dstn;
            ifc = iscell(prev_dstn_1);
            if ifc == 1
                prev_dstn_1 = cell2mat(prev_dstn_1);
            end     
            
            
            next_target_flag = 0;
            
            clear next_target
            clear next_source
            
        else
            %we have queued multiple proliferations, so we need to loop
            %through them to ensure we dont miss any
            i=next_target_row(1);
            current_target = next_target(1);
            current_source = next_source(1);
            
            prev_speed_c = prev_speed(1);
            prev_speed_1 = cell2mat(prev_speed_c);
            
            
            prev_dstn_c = prev_dstn(1);
            prev_dstn_1 = cell2mat(prev_dstn_c);
                       
            
            ks = find(spotid == current_source);
            prev_x_pos = xcoord(ks);
            prev_y_pos = ycoord(ks);
            
            prev_speed = prev_speed(2:end);
           
            prev_dstn = prev_dstn(2:end);        
            
            next_target_row = next_target_row(2:end);
            next_target = next_target(2:end);
            next_source = next_source(2:end);
        end
        
        while track_ids(i) == curr_trackid %loop through again until we reach a new main track
            
            if cellID_source(i) == current_target
                %have we found the target we are looking for?
                k = find(spotid == cellID_source(i));
                %search for daughter cell ID
                
                tempdelx = xcoord(k) - prev_x_pos;
                tempdely = prev_y_pos - ycoord(k);
                
                temp_dspl = sqrt(tempdelx.^2 + tempdely.^2);
                temp_speed = temp_dspl / 0.5 %30mins interval, turn unit into um/hour
                speed = [speed;temp_speed];
                
                delx = [delx;tempdelx];
                dely = [dely;tempdely];
     
                prev_x_pos = xcoord(k);
                prev_y_pos = ycoord(k);
                
                %then set new target
                current_target = cellID_target(i);
                current_i = i;
                i=i+1;
            else
                i=i+1;
            end
            
            if i > length(track_ids)
               break
            end
        end
        if i > length(curr_trackid)
                k = find(spotid == current_target); %last cell
                tf = frames(k);
                tempdelx = xcoord(k)-prev_x_pos;
                tempdely = prev_y_pos - ycoord(k);
                      
                delx = [delx;tempdelx];
                dely = [dely;tempdely];
                
                temp_dspl = sqrt(tempdelx.^2 + tempdely.^2);
                temp_speed = temp_dspl / 0.5 %30mins interval, turn unit into um/hour
                speed = [speed;temp_speed];
                y = [prev_speed_1;speed];    
                                
                temp_track_dstn = sqrt(delx.^2 + dely.^2) %step displacement
                sub_track_dstn = sum(temp_track_dstn); %subtrack distance
                track_dstn = sub_track_dstn + prev_dstn_1 %track distance
                                               
                track_delx = xcoord(k) - init_x_pos
                track_dely = init_y_pos - ycoord(k)
                track_dspl = sqrt(track_delx.^2 + track_dely.^2) %track displacement
                dstn_dspl_ratio = track_dstn / track_dspl  %track distance/displacement ratio
               
               
        end
    
        
        if colour(current_i)== -65281 %tip cells in subtrack
           
            ddr_m_all = [ddr_m_all;dstn_dspl_ratio];
            dstn_m_all = [dstn_m_all;track_dstn];
            dspl_m_all = [dspl_m_all;track_dspl];
            
            delx_m =[delx_m;delx];
            dely_m =[dely_m;dely];
            
            Tdelx_m_all = [Tdelx_m_all;track_delx];                     
            Tdely_m_all = [Tdely_m_all;track_dely]; 
            
            
            sm = [sm,y];
            
            speed = [];
            
            delx=[];
            dely=[];
            
        else if colour(current_i) == -16776961 %vSIV cells in subtrack
          
            ddr_b_all = [ddr_b_all;dstn_dspl_ratio];
            dstn_b_all = [dstn_b_all;track_dstn];
            dspl_b_all = [dspl_b_all;track_dspl];
            
            delx_b =[delx_b;delx];
            dely_b =[dely_b;dely];
            
            Tdelx_b_all = [Tdelx_b_all;track_delx];                     
            Tdely_b_all = [Tdely_b_all;track_dely];       
            
            sb = [sb,y];
            
            speed = [];
            
            delx=[];
            dely=[];
            
        else colour(current_i) == -26368  %branch cells in subtrack
            
            ddr_o_all = [ddr_o_all;dstn_dspl_ratio];
            dstn_o_all = [dstn_o_all;track_dstn];
            dspl_o_all = [dspl_o_all;track_dspl];
            
            delx_o =[delx_o;delx];
            dely_o =[dely_o;dely];
            
            Tdelx_o_all = [Tdelx_o_all;track_delx];                     
            Tdely_o_all = [Tdely_o_all;track_dely];
            
            so = [so,y];
            
            speed = [];
            
            delx=[];
            dely=[];

            
            end
        end
    end
end
num_cell_divisions = proliferation_counter %%count total number of proliferation

[smrow,smcol]= size(sm);
if isequal(smcol,1)                     
     sm2 = sm;
else
     sm1 = sm';  %sm contains speeds of all tracks of tip cells in one embryo
     sm1_mean = mean(sm1); %mean of all tip cell tracks; mean speed of one embryo
     sm2 = sm1_mean'; 
end
% smt = timeseries(sm2, t0:tfl);

[sbrow,sbcol]= size(sb);
if isequal(sbcol,1)                     
     sb2 = sb;
else
    sb1 = sb';  
    sb1_mean = mean(sb1); 
    sb2 = sb1_mean';
end

[sorow,socol]= size(so);
if isequal(socol,1)                     
     so2 = so;
else
     so1 = so';  
     so1_mean = mean(so1); 
     so2 = so1_mean';
end


% OUTPUT Track distance (tip, vSIV, & branch)
header1= {'dstn_tip','dstn_vSIV','dstn_branch'};
write_xls({dstn_m_all, dstn_b_all, dstn_o_all} ,Output_filename, 'Track distance', header1);

% OUTPUT Track displacement (tip, vSIV, & branch)
header2= {'dspl_tip','dspl_vSIV','dspl_branch'};
write_xls({dspl_m_all,dspl_b_all,dspl_o_all} ,Output_filename, 'Track displacement', header2);

% OUTPUT Track distance/displacement Ratio (tip, vSIV, & branch)
header3= {'ddr_tip','ddr_vSIV','ddr_branch'};
write_xls({ddr_m_all,ddr_b_all,ddr_o_all} ,Output_filename, 'Track DDratio', header3);

% OUTPUT Track deltx & dely (tip, vSIV, & branch)
header4= {'Tdelx_tip','Tdely_tip','Tdelx_vSIV','Tdely_vSIV','Tdelx_branch','Tdely_branch'};
write_xls({Tdelx_m_all,Tdely_m_all,Tdelx_b_all,Tdely_b_all,Tdelx_o_all,Tdely_o_all} ,Output_filename, 'Track Deltax_Deltay', header4);


% OUTPUT Step deltx & dely (tip, vSIV, & branch)
header5= {'delx_tip','dely_tip','delx_vSIV','dely_vSIV','delx_branch','dely_branch'};
write_xls({delx_m,dely_m,delx_b,dely_b,delx_o,dely_o} ,Output_filename, 'Step deltax_deltay', header5);


% %OUTPUT Velocity of each step (tip, vSIV, & branch)
header6= {'MeanSpeed_tip','MeanSpeed_vSIV','MeanSpeed_branch'};
write_xls({sm2,sb2,so2} ,Output_filename, 'Mean Step Velocity', header6);



%%
save ('workspace')

