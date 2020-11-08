%% In the name of Allah
close all
clear all
clc
disp('??? ???? ?????? ??????')
disp('?????? ??? ???? ')
disp(' Rezayat. et al,  Nat. Comm. 2020')
disp(' Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures ')
disp(' BCoL (IPM) 2020')
addpath(genpath('F:\Data\Reaserch\Thesis\FEF_IT project\Sample_session\Codes'))
%% Calculate Trial based data
% this code is for load data and reconfiguration of dataset
cal=0;
if cal
    patch_dat='F:\Data\Reaserch\Thesis\FEF_IT project\Sample_session\Experiment_Kopol_2018-03-18_H20_6506_L7_8000\';
    name='Experiment_Kopol_2018-03-18_H20_6506_L7_8000';
    load('F:\Data\Reaserch\Thesis\FEF_IT project\Sample_session\Eye_TF_polyfit.mat');
    Eye_TF=Eye_TF_polyfit;
    TF=Eye_TF(1).tf;
    sec= Read_trials_dms(patch_dat,TF,name);
    resp=[];
    resp(1).IT=sec.suIT;
    resp(1).FEF=sec.suFEF;
    resp(1).condition=sec.triInf(:,1);
    resp(1).bhv=sec.triInf(:,6);
    resp(1).LFP_IT=sec.LFP_IT;
    resp(1).LFP_FEF=sec.LFP_FEF;
    resp(1).saccades=sec.saccades;
    resp(1).microsaccades=sec.microsaccades;
    resp(1).eye=sec.eye;
    resp(1).name=sec.name;
    save('F:\Data\Reaserch\Thesis\FEF_IT project\Sample_session\sample_session','resp')
end
load('F:\Data\Reaserch\Thesis\FEF_IT project\Sample_session\sample_session')

%% Preprocessing
resp=Preprocessing(resp); 

%% Fig. 1b S1 Performance
Fig1b_FigS1_Performance(resp)

%% Fig. 1c, S2 a-d Firing rate 
Fig1c_Firing_rate(resp);

%% Fig. 1d S3e-f LFP Power
Fig1d_Power(resp)

%% Fig. 2 S3 PPL
Fig2_PPL(resp);

%% Fig. 3 Phase shift
Fig3_Phase_shift(resp)

%% Fig. 4 SPL
Fig4_SPL(resp);

%% Fig. 5 Object coding 
% not completed yet
% Fig5_Object_coding(resp)
