clc, 
clear all
clc

plain(1,:)=randint(1,64)
key=randint(1,64);
cipher=Encryption(plain(1,:), key);
plain(2,:)=Decryption(cipher, key)
