function varargout = SecureKeyManagement(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SecureKeyManagement_OpeningFcn, ...
                   'gui_OutputFcn',  @SecureKeyManagement_OutputFcn, ...
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

% --- Executes just before SecureKeyManagement is made visible.
function SecureKeyManagement_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for SecureKeyManagement
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% Set these in GUI initially
set(handles.popupmenu_mnode_clst, 'Enable', 'off');
set(handles.popupmenu_no_mnode, 'Enable', 'off');
set(handles.edit_x_mnode, 'Enable', 'off');
set(handles.edit_y_mnode, 'Enable', 'off');
set(handles.pushbutton_place, 'Enable', 'off');
set(handles.pushbutton_saveplots, 'Enable', 'off');
set(handles.checkbox_skm, 'Enable', 'off');
set(handles.pushbutton_run, 'Enable', 'off');

% --- Outputs from this function are returned to the command line.
function varargout = SecureKeyManagement_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Secure Key Management
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
hold on;
global place_press_count m_matrix mn_r mn_c new_entry SKM WSN_ID WSN_E
global WSN_S Ei ETX ERX Efs lambda

lambda = 1; % Should be very much less than WSNs present in given cluster

place_press_count = 0;
m_matrix = 0;
mn_r = 0; 
mn_c = 1;
new_entry = true;
SKM = false;
WSN_ID = zeros(10, 15);
WSN_E = zeros(10, 15);
WSN_S = zeros(10, 15);

% Energy parameters
Ei = 0.5;                   % Initial energy in Jouls
ETX = 16.7*0.000000001;     % Energy for transmitting single bit
ERX = 36.1*0.000000001;     % Energy for receving single bit
Efs = 1.97*0.000000001;     % Amplification co-efficient

% --- Executes on selection change in popupmenu_cluster.
function popupmenu_cluster_Callback(hObject, eventdata, handles)
global matObj sn ch_range no_of_clst neighbour_list WSN_ID WSN_E WSN_S 
global distance Ei

set(handles.popupmenu_cluster, 'Enable', 'on');
set(handles.popupmenu_mnode_clst, 'Enable', 'on');

no_of_clst = get(handles.popupmenu_cluster, 'Value');
matObj = load('data.mat');
sn = 0;
for i= 1:no_of_clst
    sn = matObj.nodes_clst(i) + sn;
end

cla % clear axes

bs_x = 50;
bs_y = 50;

% Place BS
func.placeNodes(bs_x, bs_y, '^', 1, 'k', 'g', 20);
% Naming BS
text((bs_x + 1.5), bs_y, 'BS');

% Place CH with range circle
t = 0:0.05:6.28;
ch_range = 15; % Range of communication between cluster head and neighbour nodes
for i = 1:no_of_clst
    func.placeNodes(matObj.c_head_x(i), matObj.c_head_y(i), 'mo', 1, 'k', 'r', 12);
    x1 = (matObj.c_head_x(i) + ch_range * cos(t))';
    y1 = (matObj.c_head_y(i) + ch_range * sin(t))';
    XS=[x1; x1(1)];
    YS=[y1; y1(1)];
    line(XS,YS,'color','g');
    
    % Naming CH
    str_ch = strcat('CH', num2str(i));
    text((matObj.c_head_x(i) + 1), matObj.c_head_y(i), str_ch);
end

% Wired link between BS and CH
for i = 1:no_of_clst
    line([bs_x matObj.c_head_x(i)], [bs_x matObj.c_head_y(i)],'LineWidth', 2, 'LineStyle', '-', 'color', 'r');
end

% Place SN
for i=1:sn
    func.placeNodes(matObj.x(i), matObj.y(i), 'mo', 1, 'k', 'b', 5);
end

% Finding distance and forming neighbour list
for c_h = 1:no_of_clst
    col = 1;
    for neb_node = 1:sn
        dist = sqrt((matObj.c_head_x(c_h) -  matObj.x(neb_node))^2 + (matObj.c_head_y(c_h) -  matObj.y(neb_node))^2);
        if dist <= ch_range
           distance(c_h, col) = dist;
           neighbour_list(c_h, col) = neb_node;
           col = col + 1;
        end 
    end
end

% Numbering Sensor Nodes
for c_h = 1:no_of_clst
    col = 1;
    for num = 1:size(neighbour_list, 2)
        if neighbour_list(c_h, col) ~= 0
            i = neighbour_list(c_h, col);
            str_sn = strcat('WSN', num2str(num));
            text((matObj.x(i) + 0.5), matObj.y(i), str_sn);
            col = col + 1;
        end
    end
end

% Providing IDs
% Set States and Initial energy
for c_h = 1:no_of_clst
    col = 1;
    sn_id = randi([100, 999], 1, size(neighbour_list, 2));
    for num = 1:size(neighbour_list, 2)
        if neighbour_list(c_h, col) ~= 0
            i = neighbour_list(c_h, col);
            str_id = strcat('ID=', num2str(sn_id(num)));
            text(matObj.x(i), (matObj.y(i) - 2), str_id);
            WSN_S(c_h, col) = 1; % 1 - ACTIVE
            WSN_ID(c_h, col) = sn_id(num);
            WSN_E(c_h, col) = Ei;
            col = col + 1;
        end
    end
end
assignin('base', 'WSN_ID', WSN_ID);
assignin('base', 'WSN_E', WSN_E);
assignin('base', 'WSN_S', WSN_S);
assignin('base', 'distance', distance);

% --- Executes on selection change in popupmenu_mnode_clst.
function popupmenu_mnode_clst_Callback(hObject, eventdata, handles)
set(handles.popupmenu_mnode_clst, 'Enable', 'off');
set(handles.popupmenu_no_mnode, 'Enable', 'on');
global mn_r mn_c new_entry MNODE
mn_r = mn_r + 1;
mn_c = 1;
new_entry = true;
MNODE = true;

% --- Executes on selection change in popupmenu_no_mnode.
function popupmenu_no_mnode_Callback(hObject, eventdata, handles)
set(handles.popupmenu_no_mnode, 'Enable', 'off');
set(handles.pushbutton_place, 'Enable', 'on');
set(handles.edit_x_mnode, 'Enable', 'on');
set(handles.edit_y_mnode, 'Enable', 'on');

% --- Executes on button press in pushbutton_place.
function pushbutton_place_Callback(hObject, eventdata, handles)
global matObj ch_range m_matrix mn_r mn_c place_press_count new_entry

mnode_x = str2double(get(handles.edit_x_mnode, 'String'));
mnode_y = str2double(get(handles.edit_y_mnode, 'String'));

% Find distance to check if the M-node is present in valid cluster or not
clst = get(handles.popupmenu_mnode_clst, 'Value');
d = sqrt((mnode_x - matObj.c_head_x(clst))^2 + (mnode_y - matObj.c_head_y(clst))^2);

if (isnan(mnode_x) && isnan(mnode_y))
    msgbox('Invalid Value', 'Error','error');
else
    if d <= ch_range
        place_press_count = place_press_count + 1;
        func.placeNodes(mnode_x, mnode_y, 'mo', 1, 'k', 'k', 10);
        set(handles.edit_x_mnode, 'String', '');
        set(handles.edit_y_mnode, 'String', '');
        
        % Make entry in m_matrix
        if new_entry == true      % new_entry to avoid flushing of m_matrix
            m_matrix(mn_r, mn_c) = get(handles.popupmenu_mnode_clst, 'Value');
            mn_c = mn_c + 1;
            m_matrix(mn_r, mn_c) = get(handles.popupmenu_no_mnode, 'Value');
            new_entry = false;
        end
        
        mn_c = mn_c + 1;
        m_matrix(mn_r, mn_c) = mnode_x;
        mn_c = mn_c + 1;
        m_matrix(mn_r, mn_c) = mnode_y;
    else
        msgbox('Invalid position of M-node', 'Error','error');
        set(handles.edit_x_mnode, 'String', '');
        set(handles.edit_y_mnode, 'String', '');
        set(handles.popupmenu_mnode_clst, 'Enable', 'on');
        mn_r = mn_r - 1;
    end
    
    if place_press_count == get(handles.popupmenu_no_mnode, 'Value');
        set(handles.pushbutton_place, 'Enable', 'off');
        set(handles.edit_x_mnode, 'Enable', 'off');
        set(handles.edit_y_mnode, 'Enable', 'off');
        set(handles.popupmenu_mnode_clst, 'Enable', 'on');
        set(handles.checkbox_skm, 'Enable', 'on');
        set(handles.pushbutton_run, 'Enable', 'on');
        place_press_count = 0;
    end
end
assignin('base', 'm', m_matrix)

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
set(handles.pushbutton_run, 'Enable', 'off');
set(handles.checkbox_skm, 'Enable', 'off');

global no_of_clst sn ch_range matObj MNODE WSN_S WSN_E WSN_KEY 
global m_matrix distance SKM ETX Efs  
global matrix_G matrix_K

if SKM ~= true
    q = func.primeNo(500, 1024); % 10-bit encryption
else
    q = func.primeNo(4290000000, 4294967295); % 32-bit encryption
end
assignin('base', 'q', q)

% Uses MuPAD library to calculate primitive root of given prime number
g = feval(symengine, 'numlib::primroot', q);

% Generate matrix_G
func.genarate_matrix_G(g, q)
assignin('base', 'matrix_G', matrix_G)

% BS provides keys to CH
for i = 1:no_of_clst
    line([50 matObj.c_head_x(i)], [50 matObj.c_head_y(i)],'LineWidth', 2, 'LineStyle', '-', 'color', 'g');
    pause(200/1000);
    str_key = strcat('K=', num2str(matObj.ch_key(i)));
    text(matObj.c_head_x(i), (matObj.c_head_y(i) - 2), str_key);
    line([50 matObj.c_head_x(i)], [50 matObj.c_head_y(i)],'LineWidth', 2, 'LineStyle', '-', 'color', 'r');
    pause(200/1000);
end

s_node = 1;
for rounds = 1:no_of_clst
    % matrix_D generation
    func.generate_matrix_D(10, 100, q);
    
    % matrix_A generation
    func.generate_matrix_A(rounds, q);
    
    % matrix_K generation
    func.generate_matrix_K(rounds, q);
    
    if SKM ~= true
        for s_n = 1:matObj.nodes_clst(rounds)
            func.expandingCircle(matObj.c_head_x(rounds), matObj.c_head_y(rounds), 15, 'g');
            sn_key = matrix_K(1, s_n);
            str_key = strcat('K=', num2str(sn_key));
            text((matObj.x(s_node) + 2), (matObj.y(s_node) - 2), str_key);
            WSN_KEY(rounds, s_n) = sn_key;
            s_node = s_node + 1;
        end
    end
end

if MNODE == true
    set_plot = get(handles.popupmenu_mnode_clst, 'Value');
    
    % Malicious node randomly makes some nodes IDEAL
    for i = 1:size(m_matrix, 1)
        m_clst = m_matrix(i, 1);    % Find cluster on which malicious node present
        % Make nodes IDLE
        no_idle_nodes = randi([2, 3], 1, 1); % No of IDLE nodes min - 2, max - 3
        for j = 1:no_idle_nodes
            c = randi([1, matObj.nodes_clst(set_plot)], 1, 1);
            if (WSN_S(m_clst, c) ~= 0) && (SKM == false)
                WSN_S(m_clst, c) = 2;  % 2 - IDLE
            end
        end
    end
    
    % Update axes with IDLE state
    for r = 1:no_of_clst
        for c = 1:size(WSN_S, 2)
            if WSN_S(r, c) == 1
                text(matObj.wsn_x(r, c), (matObj.wsn_y(r, c) - 4), 'ACTIVE', 'Color', 'g');
            elseif WSN_S(r, c) == 2
                func.expandingCircle(m_matrix(1, 3), m_matrix(1, 4), 15, 'k');
                text(matObj.wsn_x(r, c), (matObj.wsn_y(r, c) - 4), 'IDLE  ', 'Color', 'r');
            end
        end
    end
end
assignin('base', 'wsn_s', WSN_S);

% Attack, Residual energy, Dead nodes and plots
re_sn = sn;
col = 1;
row = 1;
rounds = 500;
dr = ch_range;   % Reference distance in meters
packets_to_ch = 0;
dropped = 0;
e_thr = 0.1;
alive_count = zeros(1, rounds);
dead_count = zeros(1, rounds);

file = false;
for round = 1:rounds
    % Done only when SKM is introduced
    if SKM == true
        if file == false
            fileID_keys = fopen('Keys.txt','w');
            fprintf(fileID_keys,'Key File: Key strength - 32 bit\n');
            fprintf(fileID_keys,'--------------------------------');
            fprintf(fileID_keys,'\n');
            file = true;
        end
        fprintf(fileID_keys, 'Round %s\n', num2str(round));
        % Freshness in keys
        for ch = 1:no_of_clst
            % matrix_D generation
            func.generate_matrix_D(10000, 100000, q);
            
            %  matrix_A generation
            func.generate_matrix_A(ch, q);
            
            % matrix_K generation
            func.generate_matrix_K(ch, q);
            
            % Save Keys in file for every iterations
            fprintf(fileID_keys, 'CH%s ', num2str(ch));
            for c = 1:matObj.nodes_clst(ch)
                str = num2str(matrix_K(1, c));
                fprintf(fileID_keys, str);
                fprintf(fileID_keys, ' ');
            end
            fprintf(fileID_keys, '\n');
        end
    end
    
    % Common task
    present_clst = 1;
    while re_sn > 0
        if matObj.wsn_x(row, col) ~= 0
            if WSN_S(row, col) == 1        % When state is ACTIVE WSN sends packets and calculates residual energy
                % Packet sent from WSN to CH
                % func.expandingCircle(matObj.wsn_x(row, col), matObj.wsn_y(row, col), 15, 'b');
                packets_to_ch = packets_to_ch + 1;
                
                if distance(row, col) <= dr && WSN_E(row, col) ~= 0 && WSN_E(row, col) > e_thr
                    WSN_E(row, col) = WSN_E(row, col) - (((ETX * 4000) + (Efs * 4000 * distance(row, col)^2))*2);
                end
            elseif WSN_S(row, col) == 2    % When state is IDLE malicious node sends packets
                % func.expandingCircle(m_matrix(1, 3), m_matrix(1, 4), 15, 'k');
                if SKM ~= true
                    dropped = dropped + 1;
                end
                
                % ACK from CH to WSN
                % func.expandingCircle(matObj.c_head_x(row), matObj.c_head_y(row), 15, 'r');
            end
            
            if present_clst == 1
                e_round1(round, col) = WSN_E(row, col);
                assignin('base', 'e_round1', e_round1);
            elseif present_clst == 2
                e_round2(round, col) = WSN_E(row, col);
                assignin('base', 'e_round2', e_round2);
            elseif present_clst == 3
                e_round3(round, col) = WSN_E(row, col);
                assignin('base', 'e_round3', e_round3);
            end
            row = row + 1; % Move to next cluster
            present_clst = present_clst + 1;
            if present_clst > set_plot
                present_clst = 1;
            end
        else
            row = row + 1; % Move to next cluster
            present_clst = present_clst + 1;
            if present_clst > set_plot
                present_clst = 1;
            end
        end
        
        if row > no_of_clst
            row = 1;
            col = col + 1;
            re_sn = re_sn - 1;
        end
        
        if no_of_clst == 1
            if col > 10
                col = 1;
                re_sn = sn;
                break
            end
        elseif col > 15
            col = 1;
            re_sn = sn;
            break
        end
    end
    
    % When completes round
    % Calculate alive node and dead node
    for i = 1:no_of_clst
        for j = 1:matObj.nodes_clst(row)
            if WSN_E(i, j) < e_thr
                dead_count(round) = dead_count(round) + 1;
            end
            alive_count(round) = sn - dead_count(round);
        end
    end
    assignin('base', 'dead_count', dead_count);
    assignin('base', 'alive_count', alive_count);
    
    % Throughput
    packets(round) = packets_to_ch;
    dropped_packets(round) = dropped;
    
    str1 = num2str((round * 100)/rounds);
    str2 = '%';
    str = strcat(str1, str2);
    set(handles.text_progress_percent, 'String', str);
end
if SKM ==true
    fclose(fileID_keys);
end
assignin('base', 'wsn_e', WSN_E);
assignin('base', 'packets', packets);
assignin('base', 'dropped_packets', dropped_packets);

ir_e1 = zeros(rounds, 1);
ir_e2 = zeros(rounds, 1);
ir_e3 = zeros(rounds, 1);
ex_e1 = zeros(rounds, 1);
ex_e2 = zeros(rounds, 1);
ex_e3 = zeros(rounds, 1);

% Count idle nodes (Exhausion Energy)
idle1 = 0;
idle2 = 0;
idle3 = 0;
switch set_plot
    case 1
        for c = 1:15
            if WSN_S(1, c) == 2
                idle1 = idle1 + 1;
            end
        end
    case 2
        for c = 1:15
            if WSN_S(1, c) == 2
                idle1 = idle1 + 1;
            end
        end
        for c = 1:15
            if WSN_S(2, c) == 2
                idle2 = idle2 + 1;
            end
        end
    case 3
        for c = 1:15
            if WSN_S(1, c) == 2
                idle1 = idle1 + 1;
            end
        end
        for c = 1:15
            if WSN_S(2, c) == 2
                idle2 = idle2 + 1;
            end
        end
        for c = 1:15
            if WSN_S(3, c) == 2
                idle3 = idle3 + 1;
            end
        end
end

switch set_plot
    case 1
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(1)
                ir_e1(i, 1) = ir_e1(i, 1) + e_round1(i, j);
                assignin('base', 'ir_e1', ir_e1);
            end
        end
        
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(1)
                if WSN_S(1, j) ~= 2
                    ex_e1(i, 1) = ex_e1(i, 1) + e_round1(i, 1);
                    assignin('base', 'ex_e1', ex_e1);
                end
            end
        end
        
        % Calulate avg
        for i = 1:rounds
            ir_e1(i, 1) = ir_e1(i, 1)/matObj.nodes_clst(1);
            ex_e1(i, 1) = ex_e1(i, 1)/(matObj.nodes_clst(1) - idle1);
            assignin('base', 'ex_e1', ex_e1);
        end
    case 2
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(1)
                ir_e1(i, 1) = ir_e1(i, 1) + e_round1(i, j);
                assignin('base', 'ir_e1', ir_e1);
            end
        end
        
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(1)
                if WSN_S(1, j) ~= 2
                    ex_e1(i, 1) = ex_e1(i, 1) + e_round1(i, 1);
                    assignin('base', 'ex_e1', ex_e1);
                end
            end
        end
        
        % Calulate avg
        for i = 1:rounds
            ir_e1(i, 1) = ir_e1(i, 1)/matObj.nodes_clst(1);
            ex_e1(i, 1) = ex_e1(i, 1)/(matObj.nodes_clst(1) - idle1);
            assignin('base', 'ex_e1', ex_e1);
        end
        
        % Calc for cluster 2
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(2)
                ir_e2(i, 1) = ir_e2(i, 1) + e_round2(i, j);
                assignin('base', 'ir_e2', ir_e2);
            end
        end
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(2)
                if WSN_S(1, j) ~= 2
                    ex_e2(i, 1) = ex_e2(i, 1) + e_round2(i, 1);
                    assignin('base', 'ex_e2', ex_e2);
                end
            end
        end
        
        % Calculate avg
        for i = 1:rounds
            ir_e2(i, 1) = ir_e2(i, 1)/matObj.nodes_clst(2);
            ex_e2(i, 1) = ex_e2(i, 1)/(matObj.nodes_clst(2) - idle2);
            assignin('base', 'ex_e2', ex_e2);
        end
    case 3
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(1)
                ir_e1(i, 1) = ir_e1(i, 1) + e_round1(i, j);
                assignin('base', 'ir_e1', ir_e1);
            end
        end
        
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(1)
                if WSN_S(1, j) ~= 2
                    ex_e1(i, 1) = ex_e1(i, 1) + e_round1(i, 1);
                    assignin('base', 'ex_e1', ex_e1);
                end
            end
        end
        
        % Calulate avg
        for i = 1:rounds
            ir_e1(i, 1) = ir_e1(i, 1)/matObj.nodes_clst(1);
            ex_e1(i, 1) = ex_e1(i, 1)/(matObj.nodes_clst(1) - idle1);
            assignin('base', 'ex_e1', ex_e1);
        end
        
        % Calc for cluster 2
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(2)
                ir_e2(i, 1) = ir_e2(i, 1) + e_round2(i, j);
                assignin('base', 'ir_e2', ir_e2);
            end
        end
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(2)
                if WSN_S(1, j) ~= 2
                    ex_e2(i, 1) = ex_e2(i, 1) + e_round2(i, 1);
                    assignin('base', 'ex_e2', ex_e2);
                end
            end
        end
        
        % Calculate avg
        for i = 1:rounds
            ir_e2(i, 1) = ir_e2(i, 1)/matObj.nodes_clst(2);
            ex_e2(i, 1) = ex_e2(i, 1)/(matObj.nodes_clst(2) - idle2);
            assignin('base', 'ex_e2', ex_e2);
        end
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(3)
                ir_e3(i, 1) = ir_e3(i, 1) + e_round3(i, j);
                assignin('base', 'ir_e3', ir_e3);
            end
        end
        
        for i = 1:rounds
            for j = 1:matObj.nodes_clst(3)
                if WSN_S(1, j) ~= 2
                    ex_e3(i, 1) = ex_e3(i, 1) + e_round3(i, 1);
                    assignin('base', 'ex_e3', ex_e3);
                end
            end
        end
        
        % Calculate avg
        for i = 1:rounds
            ir_e3(i, 1) = ir_e3(i, 1)/matObj.nodes_clst(3);
            ex_e3(i, 1) = ex_e3(i, 1)/(matObj.nodes_clst(3) - idle3);
            assignin('base', 'ex_e3', ex_e3);
        end
end

avg_ir = zeros(rounds, 1);
avg_ex = zeros(rounds, 1);

switch set_plot
    case 1
        % only 1 cluster is selected IR
        for i = 1:rounds
            if ir_e1(i, 1) >= e_thr
                index1 = i;
            end
            if ex_e1(i, 1) >= e_thr
                index2 = i;
            end
        end
        
        figure;
        subplot(1, 2, 1)
        plot(ir_e1(:,1), 'LineWidth', 2, 'Color', 'g');
        axis([0 index1 0 1])
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Energy', 'FontSize', 14)
        legend('Identity Replication')
        grid on 
        
        subplot(1, 2, 2)
        plot(ex_e1(:,1), 'LineWidth', 2, 'Color', 'r');
        axis([0 index2 0 1])
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Energy', 'FontSize', 14)
        legend('Exhaustion')
        grid on  
    case 2
        % Find avg of 2 clusters IR
        for i = 1:rounds
            avg_ir(i, 1) = (ir_e1(i, 1) + ir_e2(i, 1))/2;
            avg_ex(i, 1) = (ex_e1(i, 1) + ex_e2(i, 1))/2;
        end
        for i = 1:rounds
            if avg_ir(i, 1) >= e_thr
                index1 = i;
            end
            if avg_ex(i, 1) >= e_thr
                index2 = i;
            end
        end
        figure;
        subplot(1, 2, 1)
        plot(avg_ir(:,1), 'LineWidth', 2, 'Color', 'g');
        axis([0 index1 0 1])
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Energy', 'FontSize', 14)
        legend('Identity Replication')
        grid on
        
        subplot(1, 2, 2)
        plot(avg_ex(:,1), 'LineWidth', 2, 'Color', 'r');
        axis([0 index2 0 1])
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Energy', 'FontSize', 14)
        legend('Exhaustion')
        grid on  
    case 3
        % Find avg of 3 clusters
        for i = 1:rounds
            avg_ir(i, 1) = (ir_e1(i, 1) + ir_e2(i, 1) + ir_e3(i, 1))/3;
            avg_ex(i, 1) = (ex_e1(i, 1) + ex_e2(i, 1) + ex_e3(i, 1))/3;
        end
        for i = 1:rounds
            if avg_ir(i, 1) >= e_thr
                index1 = i;
            end
            if avg_ex(i, 1) >= e_thr
                index2 = i;
            end
        end
        figure;
        subplot(1, 2, 1)
        plot(avg_ir(:,1), 'LineWidth', 2, 'Color', 'g');
        axis([0 index1 0 1])
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Energy', 'FontSize', 14)
        legend('Identity Replication')
        grid on  
        
        subplot(1, 2, 2)
        plot(avg_ex(:,1), 'LineWidth', 2, 'Color', 'r');
        axis([0 index2 0 1])
        xlabel('Rounds', 'FontSize', 14)
        ylabel('Energy', 'FontSize', 14)
        legend('Exhaustion')
        grid on
end

% Ploting alive and dead nodes
figure;
subplot(1, 2, 1)
plot(alive_count, 'LineWidth', 2, 'Color', 'g')
xlabel('Rounds', 'FontSize', 14)
ylabel('Alive nodes', 'FontSize', 14)
grid on
subplot(1, 2, 2)
plot(dead_count, 'LineWidth', 2, 'Color', 'r')
xlabel('Rounds', 'FontSize', 14)
ylabel('Dead nodes', 'FontSize', 14)
grid on

% Measure Throughput
for i = 1:rounds
    %tp(1, i) = (packets(i)/(dropped_packets(i) + packets(i)))*100;
    tp(1, i) = packets(i) - dropped_packets(i);
end
assignin('base', 'tp', tp);
figure;
plot(tp, 'LineWidth', 2)
xlabel('Rounds', 'FontSize', 14)
ylabel('Packets received successfully', 'FontSize', 14)
legend('Throughput')
grid on
set(handles.pushbutton_saveplots, 'Enable', 'on');

% --- Executes on button press in checkbox_skm.
function checkbox_skm_Callback(hObject, eventdata, handles)
set(handles.checkbox_skm, 'Enable', 'off');
global SKM
SKM = true;

% --- Executes on button press in pushbutton_saveplots.
function pushbutton_saveplots_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_cluster_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
% --- Executes during object creation, after setting all properties.
function popupmenu_no_mnode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
function edit_x_mnode_Callback(hObject, eventdata, handles)
 
% --- Executes during object creation, after setting all properties.
function edit_x_mnode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
function edit_y_mnode_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_y_mnode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_mnode_clst_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');