% test stressPostprocessing
% Last Modified: June 1, 2016

clear all
close all
load('loadForTest')
[sigmaxx,sigmayy,sigmazz,sigmaxy,sigmaxz,sigmayz]=...
                        stressPostprocessing(uh,vh,wh,lam,mu,T,k)
                    
