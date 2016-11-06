% STEP 0: import data from txt file
% read data from txt file into cell type
filename = 'Illustration_Data.txt';
fid = fopen(filename);
data_cell = textscan(fid,'%s','delimiter','\n');
all_data_cell = data_cell{1};
empty_line_number = find(cellfun(@isempty,all_data_cell));
all_data_total_line = size(all_data_cell,1);
total_machine_cell = all_data_cell(1,:);
rows_cell = all_data_cell(2,:);
cols_cell = all_data_cell(3,:);
% % flow cost matrix cell type
FC_cell = all_data_cell((empty_line_number(1)+1):(empty_line_number(2)-1),:);
% % set up cost matrix cell tyep
SUC_cell = all_data_cell((empty_line_number(3)+1):all_data_total_line,:);
 
% convert cell type into numeric double type
% % single value
total_machines = str2double(cell2mat(total_machine_cell));
rows = str2double(cell2mat(rows_cell));
cols = str2double(cell2mat(cols_cell));
% % matrix values
FC_cell_total_row = size(FC_cell,1);
FC_cell_total_col = size(str2num(cell2mat(FC_cell(1,:))),2);
SUC_cell_total_row = size(SUC_cell,1);
SUC_cell_total_col = size(str2num(cell2mat(SUC_cell(1,:))),2);
FC = zeros(FC_cell_total_row, FC_cell_total_col);
SUC = zeros(SUC_cell_total_row,SUC_cell_total_col);
for i = 1:FC_cell_total_row
    FC(i,:) = str2num(cell2mat(FC_cell(i,:)));
end
for i = 1:SUC_cell_total_row
    SUC(i,:) = str2num(cell2mat(SUC_cell(i,:)));
end

%%%%CODE STARTS HERE

%matrix that will store the solution called floor
%matrix that will hold set up costs called final SUC
%matrix that will hold locations of machines on floor
floor = zeros(rows, cols);
finalSUC = zeros(rows,cols);
location = zeros(2, total_machines);

%find the amount of connections each machine has, 
%store in vector numarcs

    
numarcs = zeros(2,total_machines);

for i = 1:total_machines

    for j = 1:size(FC,1)
        if FC(j,1) == i
            numarcs(1,i) = numarcs(1,i) + 1;
            numarcs(2,i) = (FC(j,3) * FC(j,4)) + numarcs(2,i);
        end
        
        if FC(j,2) == i
            numarcs(1,i) = numarcs(1,i) + 1;
            numarcs(2,i) = (FC(j,3) * FC(j,4)) + numarcs(2,i);
        end
                       
    end

end


while any(numarcs(1,:) ~= 0)
%find node with most connections or, if tied, highest transportation cost
%to find priority node 

priority = 1;
cmax = 0;

for j = 1:total_machines
    if numarcs(1,j) > cmax
        priority = j;
        cmax = numarcs(1,j);    
    elseif numarcs(1,j) == cmax
        if numarcs(2,j) > numarcs(2,priority)
            priority = j;
            cmax = numarcs(1,j);
        end
    end
end
   

            


%%%place 1st machine in lowest set up space

if ~any(floor(:) == priority)  %%%%if priority is already on the board don't place it again
    
if sum(floor(:) == 0) == rows*cols %%% if priority is the first block being layed in the factory

s = 1;
minsuc = 10000;
tempsuc = 1;
tempcolfirst = 1;
temprowfirst = 1;


for i = 1:rows;
    for j = 1:cols;
        tempsuc = SUC(s,priority);
        if tempsuc <= minsuc
            minsuc= tempsuc;
            tempcolfirst = j;
            temprowfirst = i;
        end
        s= s+1;
    end    
end

floor(temprowfirst,tempcolfirst) = priority;
location(1,priority) = temprowfirst;
location(2,priority) = tempcolfirst;
finalSUC(temprowfirst,tempcolfirst) = minsuc;

else %%% for placing a priority thats not the first block,
  %%%place it in the location that has the lowest transportation costs
  %%%%%to each of the blocks that it is connected to that are already layed down
    
costadd = zeros(3, sum(floor(:) ==0));
costx = 1;
costy = 1;
count = 1;
succount = 0;
for r = 1:rows
    for c = 1:cols
        if floor(r,c) == 0
            for k = 1:size(FC,1);
                if FC(k,1) == priority
                    costadd(1,count) = FC(k,3) * FC(k,4) * abs(r - location(1,FC(k,2))) + abs(c - location(1,FC(k,2))) + costadd (1,count); 
                end
                if FC(k,2) == priority
                    costadd(1,count) = FC(k,3) * FC(k,4) * abs(i - location(1,FC(k,1))) + abs(c - location(1,FC(k,1))) + costadd(1,count); 
                end
            end
            costadd(2,count) = r;
            costadd(3,count) = c;
            count = count + 1;
        end
    end
end

minc = 10000;
a = 1;
b = 1;
for f = 1: size(costadd,2)
    if costadd(1,f) <= minc;
        minc = costadd(1,f);
        a = costadd(2,f);
        b = costadd(3,f);
    end
end
   
floor(a,b) = priority;
location(1,priority) = a;
location(2,priority) = b;
finalSUC(a,b) = SUC((a-1)*cols + b, priority);   

end
end





%%%pick which block to take next based on connection to priority block and highest
%%%transportation cost (without distance included)
z = numarcs(1,priority);

for k=1:z

transc = 0;
temptransc = 0;
next = 0;
line = 1;
    for j = 1:size(FC,1)
        if FC(j,1) == priority
            temptransc = (FC(j,3) * FC(j,4));
            if temptransc >= transc;
                transc = temptransc;
                next = FC(j,2);
                line = j;
            end
        end
        if FC(j,2) == priority
            temptransc = (FC(j,3) * FC(j,4));
            if temptransc >= transc;
                transc = temptransc;
                next = FC(j,1);
                line = j;
            end
        end
                   
    end
    
%%%%place next on board if it is not already on the board    
if any(floor(:) == next)  %%%if 
    m = FC(line,3) * FC(line,4);
    n = abs(location(1,priority) - location(1,next)) + abs(location(2,priority) - location(2,next));
    FC(line,4) = m*n;
    FC(line,3) = 0;
    numarcs(1,priority) = numarcs(1,priority) - 1;
    if numarcs(1,next) ~=0
    numarcs(1,next) = numarcs(1,next) -1;
    end
else
%%if next is not already on the board, place it in location in terms of lowest SUC + lowest tc*d    
    
s = 1;
minsuc = 10000;
tempsuc = 1;
tempcolnext = 1;
temprownext = 1;
tempdistc = 1;
temptotalc = 0;
totalc = 100000;

for i = 1:rows;
    for j = 1:cols;
        tempdistc = (abs(i - temprowfirst) + abs(j - tempcolfirst)) * transc;
        tempsuc = SUC(s,next);
        temptotalc = tempdistc + tempsuc;
        if temptotalc <= totalc && floor(i,j) ==0 
            minsuc= tempsuc;
            totalc = temptotalc;
            tempcolnext = j;
            temprownext = i;
        end
        s= s+1;
    end
end

floor(temprownext, tempcolnext) = next;
numarcs(1,priority) = numarcs(1,priority) - 1;
location(1,next) = temprownext;
location(2,next) = tempcolnext;
finalSUC(temprownext,tempcolnext) = minsuc;
 m = FC(line,3) * FC(line,4);
 n = abs(location(1,priority) - location(1,next)) + abs(location(2,priority) - location(2,next));
 FC(line,4) = m*n;
 FC(line,3) = 0;
if numarcs(1,next) ~=0
    numarcs(1,next) = numarcs(1,next) -1;
end
end
end
end

floor
final_total_cost = sum(finalSUC(:)) + sum(FC(:,4))








 
        
    
    
