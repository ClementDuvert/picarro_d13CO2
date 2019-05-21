function[Resp_dist]=respiration(sitelist,Headspace)
%
%% Function to determine the respiration rate of soil cores.
% Import data from Picarro textfiles; calculate respiration rates for whole batches of measurements.
%
% Version 2 - Nov 2017
% C. Duvert
%
% How to run this script: type respiration(sitelist,Headspace) in the command window.
%
% INPUTS
% 'sitelist' should be a cell array with character vectors
% corresponding to the Picarro file names.
% Example: sitelist={'AGG1-1-s';'GPS1-6-s';'MRS1-4-s';'RJ1-1-s';'GPP1-3-d';'MRS1-4-d'}
%
% 'Headspace' should be a single column array with the headspace volume
% for each core in dm3.
% Example: Headspace=[0.085;0.073;0.081;0.077;0.085;0.064];
%
%
% OUTPUTS
% 'Resp_dist' returns a distribution of respiration rate
% results for each core, according to Monte Carlo simulations
% (returns a 7-column array with min, q10, q25, median, q75, q90, and max).
%
% This function requires tightsubplot.m
%
%
%****************************************************************************************************************************************************
%% 1. Initialise variables and input parameters
%****************************************************************************************************************************************************

% 1A. choose the lower and upper endpoints for resp calculations
% (default 450-550ppm)
choice=questdlg('Which [CO2] lower endpoint for respiration rate calculations?','[CO2] lower endpoint','450 ppmV','500 ppmV','450 ppmV');
switch choice
    case '450 ppmV'
        start_point=450;
    case '500 ppmV'
        start_point=500;
end

choice=questdlg('Which [CO2] upper endpoint for respiration rate calculations?','[CO2] upper endpoint','500 ppmV','550 ppmV','500 ppmV');
switch choice
    case '500 ppmV'
        end_point=500;
    case '550 ppmV'
        end_point=550;
end

% 1B. choose the folder where the Picarro files are stored
S=sprintf('Enter the name of the folder where your data is stored.\nExample: /Users/XXXXX/Project1/PicarroData/Run1/');
prompt={S};
dlg_title=('Measurement batch');
num_lines=[1 105]; %size of the prompt window
def={'/Users/XXXXX/Project1/PicarroData/Run1/'};
name=inputdlg(prompt,dlg_title,num_lines,def);
folder=name{1};
clear S name

% 1C. choose the number of simulations
prompt={'Enter number of Monte Carlo simulations for each resp rate calculation'};
dlg_title='Monte Carlo';
num_lines=1;
def={'10000'};
answer_mc=inputdlg(prompt,dlg_title,num_lines,def);
mcmc=str2double(answer_mc{1});

% 1D. save plots
ch=questdlg('Would you like to save the plots as .eps files? CAUTION saving plots will overwrite any previous plot with same name.','Plot saving','YES','NO','NO');
switch ch
    case 'YES'
        pl=1;
    case 'NO'
        pl=0;
end
clear ch

% 1E. initialise some variables and constants
a=datenum(2016,1,1);
temperature=273+32; %temp in Kelvins
gas_constant=82.05; %universal gas constant in atm * mL / mol / K
molarmass=44.01; %molar mass of CO2 in g/mol
interval_resp=40; %interval for Monte Carlo simulations - in ppm
area=pi()*0.05^2; %surface area of each core in m2
tenmin = 0.006944444496185; % 10 minutes (1/24/6) for plotting purposes
Resp=cell(length(sitelist),2);
Resp(:,1)=sitelist;
Resp_dist=zeros(length(sitelist),7);

%****************************************************************************************************************************************************
%% 2. Fetch Picarro files for each measurement and create an array
%****************************************************************************************************************************************************

% 2A. fetch Picarro files
for n=1:length(sitelist)
    sitename=sitelist{n}; disp(sitename)
    filename=strcat(folder,sitename,'.dat');
    
    if exist(filename,'file') % check the nth core has been analysed in this batch
        startRow=2;
        formatSpec='%10s%28s%26f%25f%29f%13f%42f%26f%26f%26f%26f%26f%26f%27f%26f%26f%26f%25f%26f%26f%26f%26f%26f%26f%26f%26f%26f%26f%26f%26f%26f%s%[^\n\r]';
        fileID=fopen(filename,'r');
        dataArray=textscan(fileID,formatSpec,'Delimiter','\t','WhiteSpace','','HeaderLines',startRow-1,'ReturnOnError',false);
        fclose(fileID);
        
        % 2B. create reduced array with variables of use for the calculations
        % 1st column is date as datenum; 2nd column is [CO2]; 3rd column is d13C.
        Mat(:,1)=dataArray{:,3}+a;
        Mat(:,2)=dataArray{:,8};
        Mat(:,3)=dataArray{:,17};
        if length(Mat)>2000
            Mat=Mat(1:2000,:); %focus on [1:2000] interval to avoid bug, as resp rates are calculated on first few data
        else
        end
        
        
        %****************************************************************************************************************************************************
        %% 3. Calculate respiration rates and their uncertainties (Monte Carlo)
        %****************************************************************************************************************************************************
        
        % 3A. determine two intervals around the previously defined lower and upper endpoints
        Lower_start=abs(start_point-interval_resp/2-Mat(:,2)); %1st interval is centred around lower endpoint (450 ppm by default)
        Upper_start=abs(start_point+interval_resp/2-Mat(:,2));
        Lower_end=abs(end_point-interval_resp/2-Mat(:,2)); %2nd interval is centred around upper endpoint (550 ppm by default)
        Upper_end=abs(end_point+interval_resp/2-Mat(:,2));
        
        %make sure the endpoints are chosen where [CO2] is rising
        D1=smooth(diff(Mat(:,2)),50);
        D1(1:100)
        for y=2:length(Mat)
            if D1(y-1)<0.1
                Lower_start(y)=NaN; Upper_start(y)=NaN;
                Lower_end(y)=NaN; Upper_end(y)=NaN;
            end
        end
        Lower_start(1:2)=NaN; Upper_start(1:2)=NaN; Lower_end(1:2)=NaN; Upper_end(1:2)=NaN;
        
        %find the 4 corresponding indices
        [~,b2]=min(Upper_end);
        [~,b1]=min(Lower_end(1:b2));
        [~,a2]=min(Upper_start(1:b1));
        [~,a1]=min(Lower_start(1:a2)); %ensure a1 < a2 < b1 < b2
        
        
        % 3B. calculate respiration rates for a number of simulations
        %randomly select pairs of indices between the two defined intervals
        low=datasample(a1:a2,mcmc);
        up=datasample(b1:b2,mcmc);
        %calculate a respiration rate for each montecarlo loop
        ba=waitbar(0,['Calculating respiration rates based on ' num2str(mcmc) ' simulations for ' sitelist{n}]);
        for mc=1:mcmc
            waitbar(mc/mcmc,ba);
            %calculate ppm CO2 / hour (i.e. micromol CO2 / mol air / hour)
            dy=Mat(up(mc),2)-Mat(low(mc),2);
            dx=Mat(up(mc),1)*24-Mat(low(mc),1)*24;
            rate=dy/dx;
            %calculate respiration rate in mg CO2 / m2 / hour
            headvol=Headspace(n)*1000; %convert from dm3 to mL
            mol_air=1*headvol/(gas_constant*temperature); % number of moles of air in headspace
            %mol_air == [atm] * [mL] / ([atm]*[mL]/[mol]/[K] * [K]) == [mol]
            resp=mol_air*rate*molarmass/1000/area;
            %resp == [mol] * [micromol]/[mol]/[h] * [mg]/[micromol] / [m2] == [mg]/[h]/[m2]
            Resp{n,mc+1}=resp;
            clear rate headvol weight mol_air dy dx t s resp
        end
        close(ba)
        
        % 3C. extract dispersion parameters from all simulations
        T=[Resp{n,2:mcmc}];
        Resp_dist(n,1)=max(T);
        Resp_dist(n,2)=prctile(T,90);
        Resp_dist(n,3)=prctile(T,75);
        Resp_dist(n,4)=median(T);
        Resp_dist(n,5)=prctile(T,25);
        Resp_dist(n,6)=prctile(T,10);
        Resp_dist(n,7)=min(T);
        
        
        %****************************************************************************************************************************************************
        %% 4. Plotting
        %****************************************************************************************************************************************************
        
        % 4A. initialise figure
        if isnan(Resp_dist(n,:))==0  %check that there is something to plot. If not, go to next
            if n==1
                figure('pos',[300 500 350 400]) %[left bottom width height]
            else
            end
            clf
            ha=tight_subplot(2,1,[.08 0.04],[.16 .01],[.15 .015]);
            
            % 4B. raw data + result
            axes(ha(1));
            plot(Mat(:,1),Mat(:,2),'-k'); hold on
            plot(Mat(a1:b2,1),Mat(a1:b2,2),'-b','LineWidth',3);
            plot(Mat(a1:a2,1),Mat(a1:a2,2),'-r','LineWidth',2); plot(Mat(b1:b2,1),Mat(b1:b2,2),'-r','LineWidth',2);
            ylabel('CO_2 (ppmV)');
            tickstart=floor(Mat(1,1)*24*60/10)/(24*60/10); set(gca,'XTick',tickstart:tenmin:tickstart+1/(24*60/10));
            ylim([start_point-70,end_point+100]); xlim([Mat(1,1),Mat(b2,1)+1/24/12]); grid on
            datetick('x','HH:MM','keeplimits','keepticks');
            text('Position',[Mat(b2,1)+0.95/24/12,start_point-30],'String',{sitename,['{\itR} = ',num2str(Resp_dist(n,4),2),' mg CO_{2} m^{-2} h^{-1}']},'HorizontalAlignment','right','EdgeColor','white','BackgroundColor',[.9 .9 .9]);
            
            % 4C. uncertainty histogram
            axes(ha(2));
            [N,X]=hist(T,20);
            bar(X,N,'facecolor',[0.8 0.8 0.8])
            range=max(X)-min(X);
            xlabel ('{\itR} (mg CO_{2} m^{-2} h^{-1})'); ylim([0,max(N)+max(N)/10]); xlim([min(X)-range/10,max(X)+range/10]);
            text('Position',[Resp_dist(n,4)-range/6,max(N)],'String',['median = ',num2str(Resp_dist(n,4),4)],'EdgeColor','black','BackgroundColor',[.9 .9 .9]);
            text('Position',[Resp_dist(n,2)-range/8,mean(N)*1.3],'String',['{\itq_{90}} = ',num2str(Resp_dist(n,2),4)],'EdgeColor','black','BackgroundColor',[1 1 1]);
            text('Position',[Resp_dist(n,6)-range/7,mean(N)*1.3],'String',['{\itq_{10}} = ',num2str(Resp_dist(n,6),4)],'EdgeColor','black','BackgroundColor',[1 1 1]);
            
            % 4D. save the plot
            if pl==1
                saveas(gcf,strcat(folder,sitename,'_resp.eps'),'epsc')
            else
            end
            clear Diff1 Diff2 Mat Deriv resp a1 a2 b1 b2 T
        else
            clear Diff1 Diff2 Mat Deriv resp a1 a2 b1 b2 T
        end
    else
        Resp_dist(n,:)=NaN;
    end
end
end