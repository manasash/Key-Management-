function varargout = KeyMagmt_BIBD(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @KeyMagmt_BIBD_OpeningFcn, ...
    'gui_OutputFcn',  @KeyMagmt_BIBD_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

clear variables;
clc;


% --- Executes just before KeyMagmt_BIBD is made visible.
function KeyMagmt_BIBD_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
cla;
set(handles.popupmenu_num_clstr, 'Enable', 'on');
set(handles.popupmenu_num_nodes, 'Enable', 'on');
bs_x = 200;
bs_y = 200;
assignin('base','bs_x',bs_x);
assignin('base','bs_y',bs_y);
evalin('base','save myvars.mat');
% --- Outputs from this function are returned to the command line.

function varargout = KeyMagmt_BIBD_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
load ('myvars.mat', 'bs_x', 'bs_y');
cla;
% Place BS
plot(bs_x, bs_y,'^',...
    'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',20)
% Naming BS
text((bs_x + 10.5), bs_y, 'BS');
hold on
inc_reg = false;
assignin('base','inc_reg',inc_reg);
evalin('base','save myvars.mat');

% --- Executes on selection change in popupmenu_num_clstr.
function popupmenu_num_clstr_Callback(~, ~, handles)
set(handles.popupmenu_num_clstr, 'Enable', 'off');
blocks = get(handles.popupmenu_num_clstr, 'Value');
assignin('base','blocks',blocks);
evalin('base','save myvars.mat');

% --- Executes during object creation, after setting all properties.
function popupmenu_num_clstr_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_num_nodes.
function popupmenu_num_nodes_Callback(~, ~, handles)
set(handles.popupmenu_num_nodes, 'Enable', 'off');
load ('myvars.mat','blocks');
num_nodes1 = get(handles.popupmenu_num_nodes, 'Value');
num_nodes = num_nodes1 - 1;
assignin('base','num_nodes',num_nodes);
sn=blocks; % cluster head
assignin('base','sn',sn);
x=randi([0,400],1,sn); % x coordinate of a sn
y=randi([0,400],1,sn); % y coordinate of a sn
assignin('base','x',x);
assignin('base','y',y);
save ('myvars','x','y');
plot(x,y,'*b')
r22=cell(1,sn);
X22=cell(1,sn);
Y22=cell(1,sn);
for i=1:sn
    n1=num_nodes;
    radius = 30;
    angle = 2*pi*rand(n1,1);
    r1 = radius*sqrt(rand(n1,1));
    X1 = r1.*cos(angle)+ x(i);
    Y1 = r1.*sin(angle)+ y(i);
    plot(X1,Y1,'*b')
    hold on
    X22{1,i}=X1';
    Y22{1,i}=Y1';
    r22{1,i}=r1';
end
grid on
hold on
assignin('base','X1',X1);
assignin('base','Y1',Y1);
assignin('base','r1',r1);
assignin('base','radius',radius);
assignin('base','r22',r22);
assignin('base','X22',X22);
assignin('base','Y22',Y22);
evalin('base','save myvars.mat');
set(handles.pushbutton_RUN_reg, 'Enable', 'on');
set(handles.mal_clstr, 'Enable', 'on');
set(handles.inc_reg_qn, 'Enable', 'on');
set(handles.num_mal, 'Enable', 'on');

% --- Executes during object creation, after setting all properties.
function popupmenu_num_nodes_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in mal_clstr.
function mal_clstr_Callback(~, ~, handles)
mal_node_clstr = get(handles.mal_clstr,'Value');
% load myvars.mat;
load('myvars','num_nodes');
load ('myvars','blocks');
if mal_node_clstr > blocks
    msgbox('Invalid value','Error','error');
end
set(handles.mal_clstr, 'Enable', 'off');
assignin('base','mal_node_clstr',mal_node_clstr);
evalin('base','save myvars.mat');


function num_mal_Callback(hObject, ~, handles)
if (1)
    mal_node = str2double(get(hObject,'String'));
    load('myvars','radius','x','y','num_nodes','blocks','X22','Y22','r22');
    mal_node_clstr = get(handles.mal_clstr,'Value');
    if ((num_nodes + mal_node)> 45)
        msgbox('Cluster Limit is 45!','error');
    else
        ext_mal = mal_node;
        n11=ext_mal;
        radius = 30;
        angle = 2*pi*rand(n11,1);
        r10 = radius*sqrt(rand(n11,1));
        X10 = r10.*cos(angle)+ x(mal_node_clstr);
        Y10 = r10.*sin(angle)+ y(mal_node_clstr);
        plot(X10,Y10,'*r')
        hold on;
        hold on;
        for et = 1:ext_mal
            X22{1,mal_node_clstr}(num_nodes+et)=X10(et,1);
            Y22{1,mal_node_clstr}(num_nodes+et)=Y10(et,1);
            r22{1,mal_node_clstr}(num_nodes+et)=r10(et,1);
        end
        assignin('base','mal_node',mal_node);
        assignin('base','ext_mal',ext_mal);
        assignin('base','X10',X10);
        assignin('base','Y10',Y10);
        assignin('base','r10',r10);
        assignin('base','radius',radius);
        assignin('base','r22',r22);
        assignin('base','X22',X22);
        assignin('base','Y22',Y22);
        evalin('base','save myvars.mat');
    end
end

% --- Executes on button press in inc_reg_qn.
function inc_reg_qn_Callback(~, ~, ~)
load('myvars','mal_node','num_nodes');
if num_nodes+mal_node > 45
    msgbox('Cluster Limit is 45!','error');
    return;
else
    ext_mal = mal_node;
end
inc_reg = true;
assignin('base','inc_reg',inc_reg);
assignin('base','ext_mal',ext_mal);
evalin('base','save myvars.mat');

% --- Executes on button press in pushbutton_RUN_reg.
function pushbutton_RUN_reg_Callback(~, ~, handles)
set(handles.inc_reg_qn, 'Enable', 'off');
load('myvars','inc_reg','blocks','num_nodes','x','y','X22','Y22','bs_x','bs_y','ext_mal','mal_node_clstr');
load('keys_file','hashed_data_c','MN_key','CH_keys','SBD');
sn = blocks;
set(handles.pushbutton_RUN_reg, 'Enable', 'off');
if inc_reg == true
    funcs.register_nodes(blocks,num_nodes)
    linky = zeros(1,num_nodes);
    % Deployment of registered nodes
            for i=1:sn
                for ji=1:num_nodes
                        linky(1,ji)=line([x(i) X22{1,i}(ji)], [y(i) Y22{1,i}(ji)],'LineWidth', 0.25, 'LineStyle', '-', 'color', 'y');
                end
            end
        hold on
        for i = 1:sn
            funcs.draw_circle1(x(i),y(i),30,'g');
            line([bs_x x(i)], [bs_x y(i)],'LineWidth', 0.25, 'LineStyle', '-', 'color', 'r');
        end
        hold on
    hold on
        CH_keys_assigned = CH_keys(1:blocks,:);
        MN_keys_assigned = cell(1, blocks);
        hashed_data_assigned = cell(1, blocks);
        for fin = 1:blocks
            for count = 1:num_nodes
                MN_keys_assigned{1, fin} = MN_key{fin}(num_nodes);
                hashed_data_assigned{1, fin} = hashed_data_c{fin}(num_nodes);
            end
        end
        assignin('base','hashed_data_assigned',hashed_data_assigned);
        assignin('base','MN_keys_assigned',MN_keys_assigned);
        assignin('base','hashed_data_c',hashed_data_c);
        evalin('base','save myvars.mat');
        set(handles.select_cluster, 'Enable', 'on');
        num_mal = ext_mal;
        pause(2);
        rounds = 500;
        nw_size = blocks * num_nodes;
        %%%  THROUGHPUT
        TOTAL_PACKETS = num_nodes;
        total_packets = TOTAL_PACKETS;
        %%% with registration
        for round= 1:rounds
            drop1 = 0;
            dropped1(round) = drop1;
            total_packets1(round) = round * TOTAL_PACKETS;
            packets_received_1(1, round) = total_packets1(round) - drop1;
        end
       TOTAL_packets_received = packets_received_1*blocks;
        figure (1);
        title('With Node Registration Scheme');
        plot(TOTAL_packets_received, 'LineWidth', 1, 'Color', 'g');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Packets Received successfully', 'FontSize', 14)
        legend('Throughput (With node registration scheme)')
        grid on
        % % % with registration
        % Energy parameters (in Joules)
        ei = 0.5;                   % Initial energy in Joules
        etx = 16.7*0.000000001;     % Energy for transmitting single bit
        erx = 36.1*0.000000001;     % Energy for receving single bit
        efs = 1.97*0.000000001;     % Amplification co-efficient
        e_thr = 0.1;
        center = zeros(blocks);
        xcoords0 = zeros(num_nodes);
        ycoords0 = zeros(num_nodes);
        distance_x0 = cell(blocks);
        X_CRD0 = cell(blocks);
        Y_CRD0 = cell(blocks);
        for clustr = 1:blocks
            for nodeind = 1:num_nodes
                center(clustr) = x(clustr);
                xcoords0(nodeind) = X22{clustr}(nodeind);
                ycoords0(nodeind) = Y22{clustr}(nodeind);
                distance_x0{clustr}(nodeind) = sqrt((x(clustr) - xcoords0(nodeind))^2 + (y(clustr) - ycoords0(nodeind))^2);
            end
            X_CRD0{clustr} = xcoords0;
            Y_CRD0{clustr} = ycoords0;
        end
        wsn_e = 0.5 * ones(blocks,num_nodes);
        for round = 1:rounds
            dead_count = 0;
            for rowi = 1:blocks
                for coli = 1:num_nodes
                    if wsn_e(rowi, coli) > e_thr
                        wsn_e(rowi, coli) = wsn_e(rowi, coli) - (((etx * 4000) + (efs * 4000 * distance_x0{rowi}(coli)^2))*2);
                    else
                        dead_count = dead_count + 1;
                    end
                end
            end
            dead_count_round(round) = dead_count;
            alive_count_round(round) = nw_size - dead_count_round(round);
            e_rounds{1, round} = wsn_e; %%% energy of the nodes in a round
            total_energy(round) = sum(sum(e_rounds{1, round}, 2), 1); %%% total energy of the network during a round
        end
        figure (2);
        plot(total_energy, 'LineWidth', 2, 'Color', 'g');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Energy (J)', 'FontSize', 14)
        legend('With node registration scheme')
        grid on
        set(handles.num_mal_det, 'String', num_mal);
        figure (3);
        title('black hole attack');
        plot(dropped1, 'LineWidth', 1, 'Color', 'r');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Packets Dropped', 'FontSize', 14)
        legend('Black hole attack (With node registration scheme)')
        grid on
        figure (4);
        subplot(1, 2, 1)
        plot(alive_count_round, 'LineWidth', 2, 'Color', 'g');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Alive nodes', 'FontSize', 14)
        legend('With node registration scheme')
        grid on
        subplot(1, 2, 2)
        plot(dead_count_round, 'LineWidth', 2, 'Color', 'r');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Dead nodes', 'FontSize', 14)
        legend('With node registration scheme')
        grid on
end
if inc_reg ~= true
    % plotting cluster heads
        plot(x,y,'mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.17 0.66 .83],...
            'MarkerSize',7.8)
        hold on
        for i=1:sn
            str_ch = strcat('CH', num2str(i));
            text((x(i) + 1), y(i), str_ch);
        end
    % Deployment of registered nodes and establishing link between CHs and member nodes
            for i=1:blocks
                for ji=1:num_nodes
                    linky=line([x(i) X22{1,i}(ji)], [y(i) Y22{1,i}(ji)],'LineWidth', 0.25, 'LineStyle', '-', 'color', 'y');
                end
            end
            for jie=1:ext_mal
                linky(1,jie)=line([x(mal_node_clstr) X22{1,mal_node_clstr}(num_nodes+jie)], [y(mal_node_clstr) Y22{1,mal_node_clstr}(num_nodes+jie)],'LineWidth', 0.25, 'LineStyle', '-', 'color', 'y');
            end
        hold on
        for i = 1:sn
            funcs.draw_circle1(x(i),y(i),30,'g');
            line([bs_x x(i)], [bs_x y(i)],'LineWidth', 0.5, 'LineStyle', '-', 'color', 'r');
        end
        hold on
    blocks = get(handles.popupmenu_num_clstr, 'Value');
        CH_keys_assigned = CH_keys(1:blocks,:);
        MN_keys_assigned = cell(1, blocks);
        hashed_data_assigned = cell(1, blocks);
        for fin = 1:blocks
            for count = 1:num_nodes
                MN_keys_assigned{1, fin} = MN_key{fin}(num_nodes);
                hashed_data_assigned{1, fin} = hashed_data_c{fin}(num_nodes);
            end
        end
        for nod = 1:ext_mal
            MN_keys_assigned{1,mal_node_clstr}(num_nodes+nod) = MN_key{mal_node_clstr}(num_nodes+nod);
            hashed_data_assigned{1,mal_node_clstr}(num_nodes+nod) = hashed_data_c{mal_node_clstr}(num_nodes+nod);
        end
        assignin('base','hashed_data_assigned',hashed_data_assigned);
        assignin('base','MN_keys_assigned',MN_keys_assigned);
        assignin('base','hashed_data_c',hashed_data_c);
        evalin('base','save myvars.mat');
        set(handles.select_cluster, 'Enable', 'on');
        num_mal = ext_mal;
        pause(1);
        set(handles.select_cluster, 'Enable', 'on');
        num_mal = 0;
        assignin('base','num_mal',num_mal);
        evalin('base','save myvars.mat');
        pause(2);
        rounds = 500;
        TOTAL_PACKETS = num_nodes;
        total_packets = TOTAL_PACKETS;
        nw_size = blocks * num_nodes;
        %%% without registration
        packets_received_22 = zeros(1,rounds);
        drop2 = randi([ext_mal ceil((ext_mal+1)*(27/100)*num_nodes)]);
        for round= 1:rounds
            dropped2(round) = drop2*round;
            total_packets2(round) = (total_packets);
            if total_packets2 <= drop2
                packets_received_2(1, round) = 0;
            else
                packets_received_2(1, round) = (total_packets2(round) - drop2)*round;
                packets_received_22(1, round) = blocks * packets_received_2(1, round);
            end
        end
        figure (1);
        title('Without Node Registration Scheme');
        plot(packets_received_22, 'LineWidth', 1, 'Color', 'r');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Packets Received successfully', 'FontSize', 14)
        legend('Throughput (Without node registration scheme)')
        grid on
        % % % without registration
        % % Energy parameters (in Joules)
        ei = 0.5;                   % Initial energy
        etx = 16.7*0.000000001;     % Energy for transmitting single bit
        erx = 36.1*0.000000001;     % Energy for receving single bit
        efs = 1.97*0.000000001;     % Amplification co-efficient
        e_thr = 0.1;
        for clustr1 = 1:blocks
            for nodeind1 = 1:num_nodes
                center(clustr1) = x(clustr1);
                xcoords(nodeind1) = X22{clustr1}(nodeind1);
                ycoords(nodeind1) = Y22{clustr1}(nodeind1);
                distance_x{clustr1}(nodeind1) = sqrt((x(clustr1) - xcoords(nodeind1))^2 + (y(clustr1) - ycoords(nodeind1))^2);
            end
            X_CRD{clustr1} = xcoords;
            Y_CRD{clustr1} = ycoords;
        end
        wsn_e1 = 0.5 * ones(blocks,num_nodes);
        for round1 = 1:rounds
            dead_count1 = 0;
            for rowi1 = 1:blocks
%                 min1 = ext_mal;
                max = ceil((27/100)*num_nodes);
                affected1(1, rowi1) = randi([1 max]);
                wsn_e1(rowi1, 1:affected1(1, rowi1)) = 0;
                for coli1 = 1:num_nodes
                    if wsn_e1(rowi1, coli1) > e_thr
                        wsn_e1(rowi1, coli1) = wsn_e1(rowi1, coli1) - (((etx * 4000) + (efs * 4000 * distance_x{rowi1}(coli1)^2))*2);
                    else
                        dead_count1 = dead_count1 + 1;
                    end
                end
                wsn_e1(mal_node_clstr,mal_node_clstr) = 0;
            end
            dead_count_round1(round1) = dead_count1;
            alive_count_round1(round1) = nw_size - dead_count_round1(round1);
            e_rounds1{1, round1} = wsn_e1; %%% energy of the nodes in a round
            total_energy1(round1) = sum(sum(e_rounds1{1, round1}, 2), 1); %%% total energy of the network during a round
        end
        figure (2);
        plot(total_energy1, 'LineWidth', 2, 'Color', 'r');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Energy (J)', 'FontSize', 14)
        legend('Without node registration scheme')
        grid on
        figure (3);
        title('black hole attack');
        plot(dropped2, 'LineWidth', 1, 'Color', 'r');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Packets Dropped', 'FontSize', 14)
        legend('Black hole attack (Without node registration scheme)')
        grid on
        set(handles.num_mal_det, 'String', num_mal);
        figure (4);
        subplot(1, 2, 1)
        plot(alive_count_round1, 'LineWidth', 2, 'Color', 'g');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Alive nodes', 'FontSize', 14)
        legend('Without node registration scheme')
        grid on
        subplot(1, 2, 2)
        plot(dead_count_round1, 'LineWidth', 2, 'Color', 'r');
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Dead nodes', 'FontSize', 14)
        legend('Without node registration scheme')
        grid on
end


% --- Executes on selection change in select_cluster.
function select_cluster_Callback(~, ~, handles)
set(handles.select_cluster, 'Enable', 'off');
load('myvars','blocks');
% Insertion of new node
chs = get(handles.select_cluster,'Value');
assignin('base','chs',chs);
if (chs>blocks)
    msgbox('Invalid Value', 'Error','error');
end
evalin('base','save myvars.mat');
set(handles.new_node_X, 'Enable', 'on');
set(handles.new_node_Y, 'Enable', 'on');
set(handles.pushbutton_place, 'Enable', 'on');


% --- Executes during object creation, after setting all properties.
function select_cluster_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function new_node_X_Callback(~, ~, handles)
set(handles.new_node_X, 'Enable', 'off');
xq = str2double(get(handles.new_node_X,'String'));
assignin('base','xq',xq);
evalin('base','save myvars.mat');

% --- Executes during object creation, after setting all properties.
function new_node_X_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function new_node_Y_Callback(~, ~, handles)
set(handles.new_node_Y, 'Enable', 'off');
yq = str2double(get(handles.new_node_Y,'String'));
assignin('base','yq',yq);
evalin('base','save myvars.mat');

% --- Executes during object creation, after setting all properties.
function new_node_Y_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_place.
function pushbutton_place_Callback(~, ~, handles)
load('myvars','blocks','num_nodes','xq','yq','x','y','X22','Y22','bs_x','bs_y','chs','sn','inc_reg');
load('keys_file','MN_key','hashed_data_c');
set(handles.pushbutton_place, 'Enable', 'off');
% Registration process on arrival of new node
newnodeID = num_nodes+1;
if  newnodeID>45
   msgbox('No node can be added!','error');
    return;
else
    % Include new node if non-malicious
    malic = randi([1 num_nodes+1]);
    if newnodeID ~= malic
        if ~(((xq-x(chs))^2 + (yq-y(chs))^2) <= 900)
            msgbox('Invalid position!','error');
        else
            disp('valid position');
            pause(1);
            plot(xq,yq,'*k')
            if inc_reg == true
                li=line([bs_x xq], [bs_y yq],'LineWidth', 0.25, 'LineStyle', ':', 'color', 'k');
                lin=line([x(chs) xq], [y(chs) yq],'LineWidth', 0.25, 'LineStyle', ':', 'color', 'k');
                pause(2)
                delete(lin);
                delete(li);
                newnode=zeros(1,num_nodes);
                for con = 1:num_nodes
                    if newnodeID==malic
                        newnode(con) = line([xq X22{1,chs}(con)], [yq Y22{1,chs}(con)],'LineWidth', 0.25, 'LineStyle', '--', 'color', 'r');
                    else
                        newnode(con) = line([xq X22{1,chs}(con)], [yq Y22{1,chs}(con)],'LineWidth', 0.25, 'LineStyle', '-', 'color', 'g');
                        hold on
                    end
                end
                hold on
                pause(2);
                if(sn)
                    for mn = 1:num_nodes
                        delete(newnode(mn))
                    end
                end
                msgbox('Node Registration Done! New node is Non-Malicous!');
                pause(3);
                plot([x(chs) xq], [y(chs) yq],'LineWidth', 0.25, 'LineStyle', '-', 'color', 'y');
                new_node_keys_assigned = MN_key{chs}(num_nodes+1,:);
                new_node_hash_keys = hashed_data_c{num_nodes+1,:,chs};
                msgbox('Keys Assigned!');
                assignin('base','new_node_keys_assigned',new_node_keys_assigned');
                assignin('base','new_node_hash_keys',new_node_hash_keys');
                evalin('base','save myvars.mat');
            else
                plot([x(chs) xq], [y(chs) yq],'LineWidth', 0.25, 'LineStyle', '-', 'color', 'y');
                new_node_keys_assigned = MN_key{chs}(num_nodes+1,:);
                new_node_hash_keys = hashed_data_c{num_nodes+1,:,chs};
                msgbox('Keys Assigned!');
                assignin('base','new_node_keys_assigned',new_node_keys_assigned');
                assignin('base','new_node_hash_keys',new_node_hash_keys');
                evalin('base','save myvars.mat');
            end
        end
    else
        plot(xq,yq,'*k')
        msgbox('Node Registration Done! New node is Malicious!','REG');
        pause(2);
        plot(xq,yq,'or');
        return;
        end
end

% --- Executes during object creation, after setting all properties.
function num_malic_node_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(handles.num_malic_node, 'Enable', 'off');

function num_mal_det_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function num_mal_det_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function mal_clstr_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function num_mal_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end