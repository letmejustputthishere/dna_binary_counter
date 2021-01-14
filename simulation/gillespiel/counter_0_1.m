function counter_0_1()
% Simulate a two-state model of gene expression
import Gillespie.*

%% Reaction network:
%%% DNA_0
%   1. TF activation:      DNA_0 + Ara              --k_BAD_on-->                   DNA_0_BAD
%   2. TF deactivation:    DNA_0_BAD                --k_BAD_off-->                  DNA_0 + Ara

%   3. transcription:      DNA_0_BAD                --kM_BAD_flp-->                 DNA_0_BAD + mRNA_flp

%   4. translation:        mRNA_flp                 --kP_Flp-->                     mRNA_flp + Flp

%   ((5. recombination:      DNA_0 + Flp                --kR_Flp-->                     DNA_1 + Flp)) ignored for simplicity, assuming Flp is only ever inside a cell

%   6. recombination:      DNA_0_BAD + Flp          --kR_Flp-->                     DNA_1_BAD + Flp

%   7. mRNA_flp decay:     mRNA_flp                 --gM_flp-->                     0
%   8. Flp decay:          Flp                      --gP_Flp-->                     0
%   9. Ara decay:          Ara                      --gP_Ara-->                     0

%% Rate constants
p.k_BAD_on = 0.0002; % TF activation: sec^-1 -> calculated by hand
p.k_BAD_off = 0.01; %TF deactivation: sec^-1 -> given by kobi
p.kM_BAD_flp = 0.000231; %transcription:  sec^-1 -> calculated by hand ( steady state concentration: 0.1 ; half-life: 5*60 sec)
p.kP_Flp = 0.00577; % translation: sec^-1 -> calculated by hand ( steady state concentration: 10 ; half-life: 20*60 sec (cell divison rate))
% p.kR_Flp = 4.81*10^-5; % recombination: sec^1 -> calculated by hand ( from paper suggestiong 4h half-time for the process)
p.kR_Flp = 2.8*10^-3; % recombination: sec^1 -> derived from paper stating cleavage rate constant for Flp
p.gM_flp =0.00231; % mRNA_flp decay: sec^-1 -> calculated by hand (derived from k_syn)
p.gP_Flp  = 0.011; % Flp decay: sec^-1 -> calculated by hand (derived from Paper stating GFPssrA tagged half time of 60s)
p.gP_Ara = 0.0003; % from Paper, constant consumption rate of arabinose
% p.gP_Ara = 0.1201; % from Paper, exponential degradation 
%% Initial state
tspan = [0, 60*60*8]; %seconds (8 hour Ara pulse described in paper)
x0    = [100, 8000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];     %DNA_0, Ara, DNA_0_BAD, mRNA_flp, Flp, DNA_1, DNA_1_BAD, mRNA_gfp, mRNA_iptg, GFP, IPTG, DNA_1_BAD_A1lacO, DNA_0_BAD_A1lacO

%% Specify reaction network
pfun = @propensities_2state;
                %DNA_0, Ara, DNA_0_BAD, mRNA_flp, Flp, DNA_1, DNA_1_BAD, mRNA_gfp, mRNA_iptg, GFP, IPTG, DNA_1_BAD_A1lacO, DNA_0_BAD_A1lacO
stoich_matrix = [-1     -1      1          0       0     0        0         0          0       0     0          0               0            %DNA_0 + Ara              --k_BAD_on-->                   DNA_0_BAD
                  1      1     -1          0       0     0        0         0          0       0     0          0               0            %DNA_0_BAD                --k_BAD_off-->                  DNA_0 + Ara
                  0      0      0          1       0     0        0         0          0       0     0          0               0            %DNA_0_BAD                --kM_BAD_flp-->                 DNA_0_BAD + mRNA_flp
                  0      0      0          0       1     0        0         0          0       0     0          0               0            %mRNA_flp                 --kP_Flp-->                     mRNA_flp + Flp
                  0      0     -1          0       0     0        1         0          0       0     0          0               0            %DNA_0_BAD + 4Flp         --kR_Flp-->                     DNA_1_BAD + 4Flp
                  0      0      0         -1       0     0        0         0          0       0     0          0               0            %mRNA_flp                 --gM_flp-->                     0
                  0      0      0          0      -1     0        0         0          0       0     0          0               0            %Flp                      --gP_Flp-->                     0
                  0     -1      0          0       0     0        0         0          0       0     0          0               0];          %Ara                      --gP_Ara-->                     0


                

%% Run simulation
[t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
%[t,x] = firstReactionMethod(stoich_matrix, pfun, tspan, x0, p);

%% Plot time course
figure();

% REMOVE ARA COLUMN FROM VECTOR

% CAREFUL, INDEXES SHIFT WITH THIS OPERATION!
x(:,2) = [];

% REMOVE DNA_1 COLUMN FROM VECTOR ( NOW AT COLUMN 5 INSTEAD OF 6)
x(:,5) = [];


stairs(t,x(:,1:5)); set(gca,'XLim',tspan);

% set log scale 
set(gca,'XScale','log');

xlabel('time (s)');
ylabel('molecules');

% CAUTION: DNA_1 and Ara ARE MISSING FROM THIS GRAPH
legend('DNA_0', 'DNA_0_BAD', 'mRNA_flp', 'Flp', 'DNA_1_BAD');

end


function a = propensities_2state(x, p)
% Return reaction propensities given current state x
DNA_0    = x(1);
Ara = x(2);
DNA_0_BAD = x(3);
mRNA_flp = x(4);
Flp = x(5);
DNA_1 = x(6);
DNA_1_BAD =x(7);
mRNA_gfp = x(8);
mRNA_iptg = x(9);
GFP = x(10);
IPTG = x(11);
DNA_1_BAD_A1lacO = x(12);
DNA_0_BAD_A1lacO = x(13);

R_tot = p.k_BAD_on*DNA_0*Ara + ...
    p.k_BAD_off*DNA_0_BAD + ...
    p.kM_BAD_flp*DNA_0_BAD + ...
    p.kP_Flp*mRNA_flp + ...
    p.kR_Flp*DNA_0_BAD*Flp + ...
    p.gM_flp*mRNA_flp + ...
    p.gP_Flp*Flp + ...
    p.gP_Ara*Ara;

a = [p.k_BAD_on*DNA_0*Ara/R_tot;            %activation promoter
     p.k_BAD_off*DNA_0_BAD;       %deactivation promoter
     p.kM_BAD_flp*DNA_0_BAD;       %transcription
     p.kP_Flp*mRNA_flp;   %translation
     p.kR_Flp*DNA_0_BAD*Flp/R_tot;   %recombination
     p.gM_flp*mRNA_flp;   %protein decay
     p.gP_Flp*Flp; %protein decay 
     p.gP_Ara];   %protein decay !!! CONSTANT LIKE IN PAPER −cAra !!!
%      p.gP_Ara*Ara];   %protein decay !!! exponential −dAra ⋅[ara] !!!
 
end
