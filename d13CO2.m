function[Int_dist]=d13CO2(sitelist)
%
%% Function to determine the d13C of respired CO2 from soil cores using Keeling plots.
% Import data from Picarro textfiles; create Keeling plots for whole batches of measurements;
% calculate uncertainties for each Keeling intercept.
%
% Version 2.1 - Feb 2018
% C. Duvert - J. Wynn
%
%
% This version does not use a different calibration for each batch,
% but rather uses one "master" calibration for all batches
% (to be determined at initial step). The true standard d13C values can be
% modified in section [1E] of this script.
%
% How to run the script: type d13CO2(sitelist) in the command window. Pop-up
% windows will then guide you through the simulations.
%
% INPUTS
% 'sitelist' should be a cell array with character vectors
% corresponding to the Picarro file names.
% Example: sitelist={'AGG1-1-s';'GPS1-6-s';'MRS1-4-s';'RJ1-1-s';'GPP1-3-d';'MRS1-4-d'}
%
% OUTPUTS
% 'Int_dist' returns a distribution of Keeling intercept
% results for each core, according to Monte Carlo simulations
% (returns a 7-column array with min, q10, q25, median, q75, q90, and max).
%
% This function requires tightsubplot.m
%
%
%% 1. Initialise variables and input parameters
%***********************************************

% 1A. choose the upper endpoint
choice=questdlg('Which [CO2] upper endpoint?','[CO2] endpoint','700 ppmV','800 ppmV','700 ppmV');
switch choice
    case '700 ppmV'
        endpoint=700;
    case '800 ppmV'
        endpoint=800;
end

% 1B. choose the smoothing interval
prompt={'Enter smoothing interval for the calculation of derivatives'};
dlg_title='Smoothing';
num_lines=1;
def={'100'};
answer=inputdlg(prompt,dlg_title,num_lines,def);
smoothing=str2double(answer{1});

% 1C. choose the number of simulations
prompt={'Enter number of Monte Carlo simulations for each Keeling calculation'};
dlg_title='Monte Carlo';
num_lines=1;
def={'10000'};
answer_mc=inputdlg(prompt,dlg_title,num_lines,def);
mcmc=str2double(answer_mc{1});

% 1D. enter the folder where the Picarro files are stored
S=sprintf('Enter the name of the folder where your data is stored.\nExample: /Users/XXXXX/Project1/PicarroData/Run1/');
prompt={S};
dlg_title=('Measurement batch');
num_lines=[1 105]; %size of the prompt window
def={'/Users/XXXXX/Project1/PicarroData/Run1/'};
name=inputdlg(prompt,dlg_title,num_lines,def);
folder=name{1};
clear S name

% 1E. initialise variables and constants
Int_sim=cell(length(sitelist),mcmc+1); %create array to contain all simulations of Keeling intercepts
Int_sim(:,1)=sitelist;
Int_dist=zeros(length(sitelist),7); %create array to contain a summary of distributions
a=datenum(2016,1,1); %initialise date
interval=50; %interval for Monte Carlo simulations - in ppm
Std_gas=-38.53; %actual std values; instrument-specific and to be modified accordingly.
Std_NaHCO3=-6.02;

% 1F. set the calibration
prompt={'Enter the average measured value for lower standard (gas)'};
dlg_title='Calibration - lower standard';
num_lines=1;
def={'-35.8494'};
answer2=inputdlg(prompt,dlg_title,num_lines,def);
lowerstd=str2double(answer2{1});

prompt={'Enter the average measured value for upper standard (NaHCO3)'};
dlg_title='Calibration - upper standard';
num_lines=1;
def={'-9.0108'};
answer3=inputdlg(prompt,dlg_title,num_lines,def);
upperstd=str2double(answer3{1});

alpha=(Std_NaHCO3-Std_gas)/(upperstd-lowerstd);
beta=Std_NaHCO3-alpha*upperstd;

% 1G. save plots (optional)
choice3=questdlg('Would you like to save the plots as .eps files? CAUTION saving plots will overwrite any previous plot with same name.','Plot saving','YES','NO','NO');
switch choice3
    case 'YES'
        pl=1;
    case 'NO'
        pl=0;
end


%% 2. Fetch Picarro files for each measurement, create a reduced, corrected array and trim off its ends
%*******************************************************************************************************

% 2A. fetch Picarro files
for n=1:length(sitelist)
    sitename=sitelist{n};
    disp(sitename)
    filename=strcat(folder,sitename,'.dat');
    
    if exist(filename,'file') % check if the nth core has been analysed in this batch 
        startRow=2;
        formatSpec='%10s%28s%26f%25f%29f%13f%42f%26f%26f%26f%26f%26f%26f%27f%26f%26f%26f%25f%26f%26f%26f%26f%26f%26f%26f%26f%26f%26f%26f%26f%26f%s%[^\n\r]';
        fileID=fopen(filename,'r');
        dataArray=textscan(fileID,formatSpec,'Delimiter','\t','WhiteSpace','','HeaderLines',startRow-1,'ReturnOnError',false);
        fclose(fileID);
    
    % 2B. create reduced array with variables of use for the calculations
    % 1st column is date as datenum; 2nd column is [CO2]; 3rd column is d13C.
    K(:,1)=dataArray{:,3}+a;
    K(:,2)=dataArray{:,8};
    K(:,3)=dataArray{:,17};
    
    % 2C. correct values according to calibration
    K(:,3)=alpha*K(:,3)+beta; %3rd column is now *corrected* d13C.
    clearvars filename startRow formatSpec fileID dataArray ans;
    
    % 2D. trim off the dataset: lower endpoint (inflection point after initial bump)
    d1=diff(K(:,2)); %calculate first derivative
    av_d1=smooth(d1,smoothing); %apply low-pass filter to get rid of the noise
    d2=diff(av_d1); %calculate second derivative
    av_d2=smooth(d2,smoothing); %low-pass filter again
    [~,s]=(min(abs(av_d2(50:200)))); %find the inflection point, i.e. where the 2nd derivative is zero (50:200 is used to focus on the interval of interest)
    K_lt=K(s:end,:); %trim the dataset accordingly
    
    % 2E. trim off the dataset: upper endpoint
    Diff=abs(endpoint-K_lt(:,2));
    for y=2:length(K_lt) %loop to make sure the endpoint is chosen where [CO2] is rising
        if av_d1(y-1)<0
          if y>50
              Diff(y-50:y)=NaN; %need to assign NaNs for 50 anterior values as well due to smoothing effects
          else
              Diff(y)=NaN;
          end
        end
    end
    [~,t]=(min(Diff));
    K_t=K_lt(1:t,:);
    clear Diff
    
        
    %% 3. Quantify the uncertainty on Keeling intercepts using a number of randomly selected lower and upper endpoints
    %******************************************************************************************************************
   
   % 3A. define the intervals where to pick random endpoints
   v=t+s-1;
   if K(v,2)-K(s,2)>3*interval %rule to check that total data range is > 3 intervals. If not, reduce interval
       [~,o]=min(abs((K(s,2)+interval-K(1:v,2))));
       [~,u]=min(abs((K(1:v,2)-(K(v,2)-interval))));
   else
       interval2=floor((K(v,2)-K(s,2))/3); %reduced interval should be a third of the total data range
       [~,o]=min(abs((K(s,2)+interval2-K(1:v,2))));
       [~,u]=min(abs((K(1:v,2)-(K(v,2)-interval2))));
   end
   clear interval2
   
   % 3B. check that intervals are not negative. If so, warning message and returns NaNs
   if o-s<2 || v-u<2 
       Int_dist(n,:)=NaN;
       warningMessage=sprintf('Warning: there was an issue with the calculation of lower and/or upper endpoints for this core:\n%s', sitename);
       uiwait(msgbox(warningMessage));
       clear K fits f int y s t v o u b c T
   else

   % 3C. randomly select pairs of lower and upper endpoints
   [~,b]=datasample(K(s:o,2),mcmc);
   [~,c]=datasample(K(u:v,2),mcmc);
   low=s+b; up=u+c-1;
   
   % 3D. calculate the fit for each simulation
   ba=waitbar(0,['Calculating uncertainty based on ' num2str(mcmc) ' simulations for ' sitelist{n}]);
   for m=1:mcmc
       waitbar(m/mcmc,ba);
       K_mc=K(low(m):up(m),:);
       Rec=1./K_mc(:,2);
       [f,~]=fit(Rec,K_mc(:,3),fittype('poly1'));
       fits=coeffvalues(f); int=fits(2);
       Int_sim{n,m+1}=int;
       clear f int
   end
   close(ba)
   
   % 3E. extract distribution of error for all simulations (will be used for whisker plots)
   T=[Int_sim{n,2:mcmc}];
   Int_dist(n,1)=max(T);
   Int_dist(n,2)=prctile(T,90);
   Int_dist(n,3)=prctile(T,75);
   Int_dist(n,4)=median(T);
   Int_dist(n,5)=prctile(T,25);
   Int_dist(n,6)=prctile(T,10);
   Int_dist(n,7)=min(T);

   
    %% 4. Plotting
    %**************
    
    % 4A. initialise figure
    if n==1
        figure('pos',[10 100 1000 900]) %[left bottom width height]
    else
    end
    clf
    ha=tight_subplot(3,2,[.03 0.04],[.08 .01],[.06 .015]);
    
    % 4B. [CO2] data with overlay of trimmed data
    axes(ha(1));
    plot(K(:,1),K(:,2),'-k');
    hold on
    plot(K_t(:,1),K_t(:,2),'-b','LineWidth',3);
    plot(K(s:o,1),K(s:o,2),'-r','LineWidth',2); plot(K(u:v,1),K(u:v,2),'-r','LineWidth',2);
    datetick('x','HH:MM','keeplimits');
    ylabel ('CO_2 (ppmV)'); ylim([300,900]); xlim([K(1,1),K_t(end,1)+1/24/12]); grid on
    legend('all data','trimmed data used for Keeling plot','intervals used for uncertainty quantification','Location','northwest')
    
    axes(ha(2));
    delete(ha(2))
    
    % 4C. d13C trimmed data
    axes(ha(3));
    plot(K(:,1),K(:,3),'-k');
    hold on
    plot(K_t(:,1),K_t(:,3),'-b');
    datetick('x','HH:MM','keeplimits');
    y=K_t(:,3);
    ylabel ('\delta^{13}C (‰)'); ylim([min(y)-1,max(y)+1]); xlim([K(1,1),K_t(end,1)+1/24/12]); grid on
    
    axes(ha(4));
    delete(ha(4))
    
    % 4D. Keeling plot (based on initial interval)
    axes(ha(5));
    R_CO2=1./K(:,2); % calculate reciprocal to draw line based on initial interval 
    R_CO2_t=1./K_t(:,2);
    [f,~]=fit(R_CO2_t,K_t(:,3),fittype('poly1'));
    fits=coeffvalues(f);
    xA=min(R_CO2_t); yA=fits(1)*xA+fits(2); %determine two points to draw the line
    xB=max(R_CO2_t); yB=fits(1)*xB+fits(2);
    plot(R_CO2,K(:,3),'.k','color',[.6 .6 .6])
    hold on
    plot(R_CO2_t,K_t(:,3),'ok','MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',4)
    plot([xA xB],[yA yB],'-r','LineWidth',2);
    hold on
    xlabel('1/CO_2 (ppmV)'); ylabel('\delta^{13}C (‰)'); xlim([0.00115,0.00245]); ylim([min(y)-1,max(y)+1]);
    text('Position',[0.00118,max(y)-0.1],'String',{sitename,['Intercept = ',num2str(Int_dist(n,4),4),'‰']},'color',[1 1 1],'EdgeColor','black','BackgroundColor',[.3 .3 .3]);
    grid on
    
    % 4E. uncertainty histogram
    axes(ha(6))
    [N,X]=hist(T,20);
    bar(X,N,'facecolor',[0.8 0.8 0.8])
    range=max(X)-min(X);
    xlabel ('\delta^{13}C (‰)'); ylim([0,max(N)+max(N)/10]); xlim([min(X)-range/10,max(X)+range/10]);
    text('Position',[Int_dist(n,4)-range/6,max(N)],'String',['median = ',num2str(Int_dist(n,4),4),'‰'],'EdgeColor','black','BackgroundColor',[.9 .9 .9]);
    text('Position',[Int_dist(n,2)-range/8,mean(N)*1.7],'String',['{\itq_{90}} = ',num2str(Int_dist(n,2),4),'‰'],'EdgeColor','black','BackgroundColor',[1 1 1]);
    text('Position',[Int_dist(n,6)-range/7,mean(N)*1.7],'String',['{\itq_{10}} = ',num2str(Int_dist(n,6),4),'‰'],'EdgeColor','black','BackgroundColor',[1 1 1]);

    % 4F. save plot
    if pl==1
        saveas(gcf,strcat(folder,sitename,'.eps'),'epsc')
    else
    end
    clear K fits f int y xA xB yA yB s t v o u b c T
    
    %proceed to next plot
         ch=questdlg('Proceed to next graph?','Check','YES','NO','YES');
         switch ch
             case 'YES'
                 continue
             case 'NO'
                 pause (5)
     
        end
   end
   
    else
        Int_dist(n,:)=NaN;
        warningMessage=sprintf('Warning: this core has not been analysed:\n%s', sitename);
        uiwait(msgbox(warningMessage));
    end
end


% 4G. save intercept values as text files
txtfile=fopen('Intercepts_distrib.txt','w');
fprintf(txtfile,'%s %6s %6s %6s %6s %6s %6s %6s\n','CORE','max','q90','q75','median','q25','q10','min');
for e=1:length(sitelist)
    fprintf(txtfile,'%s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',sitelist{e},Int_dist(e,:));
end
fclose(txtfile);

return

end

