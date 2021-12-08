clc
clear all
load MPC_cos.mat

v11=sqrt(mean((yrliftstack(1,:) - XX_koop(1,1:500)).^2))
v12=sqrt(mean((yrliftstack(1,:) - Y_new_p(1,:)).^2))
v21=100*norm(yrliftstack(1,:) - XX_koop(1,1:500))/norm(yrliftstack(1,:))
v22 =100*norm(yrliftstack(1,:) - Y_new_p(1,:))/norm(yrliftstack(1,:))