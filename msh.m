% tic;
%function [MN_key, hashed_data_c] = key_generation

CH_keys = gen_bibd;

[M,N] = size(CH_keys);
% Preserve the row indices
rowIndex = repmat((1:M)',[1 N]);

% Get randomized column indices by sorting a second random array
[~,randomizedColIndex] = sort(rand(M,N),2);

% Need to use linear indexing to create B
newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
B = CH_keys(newLinearIndex);

%random sufflle of all the rows
%random_B = B(randperm(size(B, 1)), :);

for i=1:M
    C(i,1:2) = B(i,1:2);
    hmac_key(i,:) = C(i,1).*C(i,2);
    other_key(i,:) = B(i,3);
end

HMAC_key = num2str(hmac_key);
authentication_key = sym(hmac_key); % 256 bit HMAC key



for ij=1:M
    ac= C(ij,1).*other_key(ij);
    bc= C(ij,2).*other_key(ij);
    
    prod = ac.*bc;
    strprod =string(dec2bin(prod));
    
    no_of_bits=ceil(strlength(strprod)/2);
    extract256=extractBetween(strprod,1,no_of_bits);
    
    str1=sprintf('text2int("%s",2)',extract256);
    
    s_key(ij,:) = evalin(symengine,str1);
    
end

secret_key = s_key; % secret_key of 256 bit length

%hashed_data_c = cell(1,M);
SBD = cell(1,12);
MN_key = cell(1,M);
for chn = 1:M
    %keys for member nodes
    [IM,SB] = gen_sbibd;
    
    % all the member_node key of first cluster head
    %sbibd_key = secret_key(1).*XX;
    
    % displaying sbibd keys
    
    for si = 1:45
        for sj =1:12
            sbibd_k(si,sj) = sym(secret_key(chn))*sym(SB(si,sj));
        end
    end
    
    MN_key{1,chn} = sym(sbibd_k); %  382 bit member node key
    SBD{1,chn} = SB;
    
    
    % DATA  AUTHENTICATION USING HMAC SHA-256 ALGORITHM
    
%     hashed_data = cell(45,12);
%     for ir=1:45
%         for ic=1:12
%             hashed_data{ir,ic} = HMAC(HMAC_key(ic,1),'HI','SHA-256');
%         end
%     end
%     
%     hashed_data_c{1,chn} = hashed_data;
    
end

% DATA  AUTHENTICATION USING HMAC SHA-256 ALGORITHM

%Assuminng that same data that is 'HI' is processed, Data authetication
%for MN in each of the CH is shown 

hashed_data = cell(45,12,12);
 for cl = 1:12
       hashed_data(:,:,cl) = cellstr(HMAC(HMAC_key(cl,1),'HI','SHA-256'));
 end
 hashed_data_c = hashed_data;
 
 assignin('base','hashed_data_c',hashed_data_c);
 assignin('base','MN_key',MN_key);
 assignin('base','CH_keys',CH_keys);
 assignin('base','SBD',SBD);
 evalin('base','save keys_file.mat');
%  load keys_file.mat;
%  
%  save('keys_file','hashed_data_c','MN_key','CH_keys','SBD');
 
 %end
 % toc;