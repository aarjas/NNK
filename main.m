%%
load('C:/Users/aarjas/OneDrive - Oulun yliopisto/Tracking paper - reresubmission/submission/Codes+data/trackingdata.mat')
load('C:/Users/aarjas/OneDrive - Oulun yliopisto/Tracking paper - reresubmission/submission/Codes+data/net3.mat')
%%
track=nnk(images,RFdata,net3,true)