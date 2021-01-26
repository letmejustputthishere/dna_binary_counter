function [x] = counter_1_0_no_pulse(x0)
% Simulate a two bit DNA based counter
import Gillespie.*
import Figures.*

%% Reaction network:
%%% DNA_0
%   1. TF activation:      DNA_0 + Ara              --k_BAD_on-->                   DNA_0_BAD
%   2. TF deactivation:    DNA_0_BAD                --k_BAD_off-->                  DNA_0 + Ara

%   ((3. transcription:      DNA_0_BAD                --kM_BAD_flp-->                 DNA_0_BAD + mRNA_flp))

%   4. translation:        mRNA_flp                 --kP_Flp-->                     mRNA_flp + Flp

%   ((5. recombination:      DNA_0 + Flp                --kR_Flp-->                     DNA_1 + Flp)) ignored for simplicity, assuming Flp is only ever inside a cell

%   6. recombination:      DNA_0_BAD + Flp          --kR_Flp-->                     DNA_1_BAD + Flp

%   7. mRNA_flp decay:     mRNA_flp                 --gM_flp-->                     0
%   8. Flp decay:          Flp                      --gP_Flp-->                     0
%   9. Ara decay:          Ara                      --gP_Ara-->                     0

%   10. recombination:     DNA_1_BAD + Flp          --kR_Flp-->                     DNA_0_BAD + Flp

%   11. TF activation:     DNA_1 + Ara              --k_BAD_on-->                   DNA_1_BAD
%   12. TF deactivation:   DNA_1_BAD                --k_BAD_off-->                  DNA_1 + Ara

%%% DNA_1
%   13. transcription:     DNA_1_BAD                --kM_BAD_gfp_T7p-->             DNA_1_BAD + mRNA_gfp + mRNA_T7p

%   14. translation:       mRNA_gfp                 --kP_GFP-->                     mRNA_gfp + GFP
%   15. translation:       mRNA_T7p                 --kP_T7p-->                     mRNA_T7p + T7p

%   16. T7p transcription: DNA_1 + T7p              --kM_T7p_flp-->                 DNA_1 + T7p + mRNA_flp
%   17. T7p transcription: DNA_1_BAD + T7p          --kM_T7p_flp-->                 DNA_1_BAD + T7p + mRNA_flp 
%   ((18. T7p transcription: DNA_0 + T7p              --kM_T7p_flp-->                 DNA_0 + T7p + mRNA_flp)) ignored for simplicity, while the recombinase is working no more T7p can be produced. T7p has a 50 min half life so 
%   ((19. T7p transcription: DNA_0_BAD + T7p          --kM_T7p_flp-->                 DNA_0_BAD + T7p + mRNA_flp)) no T7p should be available if the DNA flipped back to state 0


%   20. recombination:     DNA_1 + Flp              --kR_Flp-->                     DNA_0 + Flp
%   21. recombination:     DNA_0 + Flp              --kR_Flp-->                     DNA_1 + Flp


%   22. mRNA_gfp decay:    mRNA_gfp                 --gM_gfp-->                     0
%   23. mRNA_T7p decay:    mRNA_T7p                 --gM_T7p-->                     0
%   24. GFP decay:         GFP                      --gP_GFP-->                     0
%   25. T7p decay:         T7p                      --gP_T7p-->                     0


%% Rate constants
p.k_BAD_on = 0.0002; % TF activation: sec^-1 -> calculated by hand
p.k_BAD_off = 0.01; %TF deactivation: sec^-1 -> given by kobi
p.kM_BAD_flp = 0.00231; %transcription:  sec^-1 -> calculated by hand ( steady state concentration: 1 ; half-life: 5*60 sec)
p.kP_Flp = 0.0577; % translation: sec^-1 -> calculated by hand ( steady state concentration: 100 ; half-life: 20*60 sec (cell divison rate))
p.kR_Flp = 4.81*10^-5; % recombination: sec^1 -> calculated by hand ( from paper suggestiong 4h half-time for the process)
% p.kR_Flp = 2.8*10^-3; % recombination: sec^1 -> derived from paper stating cleavage rate constant for Flp
p.gM_flp =0.00231; % mRNA_flp decay: sec^-1 -> calculated by hand (derived from half life )
p.gP_Flp  = 0.011; % Flp decay: sec^-1 -> calculated by hand (derived from Paper stating GFPssrA tagged half time of 60s)
%p.gP_Ara = 0.0003; % from Paper, constant consumption rate of arabinose
p.gP_Ara = 0.1201; % from Paper, exponential degradation 
p.kM_BAD_gfp_T7p = 0.00231; %transcription:  sec^-1 -> calculated by hand ( steady state concentration: 1 ; half-life: 5*60 sec) 
p.kP_GFP = 0.0577; % translation: sec^-1 -> calculated by hand ( steady state concentration: 100 ; half-life: 20*60 sec (cell divison rate))
p.kP_T7p = 0.0577; % translation: sec^-1 -> calculated by hand ( steady state concentration: 100 ; half-life: 20*60 sec (cell divison rate))
p.kM_T7p_flp = 5*0.00231; %transcription:  sec^-1 -> calculated by hand ( steady state concentration: 1 ; half-life: 5*60 sec)  (5 times higher than normal RNA polymerase)
p.gM_gfp = 0.00231; % mRNA decay: sec^-1 -> calculated by hand (derived from half life )
p.gM_T7p = 0.00231; % mRNA decay: sec^-1 -> calculated by hand (derived from half life )
p.gP_GFP = 0.000577; % calcualted by hand (dervied from 20 min half-life)
%p.gP_T7p = 0.00023; % derived from paper 
p.gP_T7p = 0.011; % T7p decay: sec^-1 -> calculated by hand (derived from Paper stating GFPssrA tagged half time of 60s)


%% Initial state
tspan = [0, 60*60*100]; %seconds (8 hour Ara pulse described in paper)


%% Specify reaction network
pfun = @propensities_2state;
                %DNA_0, Ara, DNA_0_BAD, mRNA_flp, Flp, DNA_1, DNA_1_BAD, mRNA_gfp, mRNA_T7p, GFP,   T7p 
stoich_matrix = [-1     -1      1          0       0     0        0         0          0       0     0            %DNA_0 + Ara              --k_BAD_on-->                   DNA_0_BAD
                  1      1     -1          0       0     0        0         0          0       0     0            %DNA_0_BAD                --k_BAD_off-->                  DNA_0 + Ara
                  0      0      0          0       1     0        0         0          0       0     0            %mRNA_flp                 --kP_Flp-->                     mRNA_flp + Flp
                  0      0     -1          0       0     0        1         0          0       0     0            %DNA_0_BAD + 4Flp         --kR_Flp-->                     DNA_1_BAD + 4Flp
                  0      0      0         -1       0     0        0         0          0       0     0            %mRNA_flp                 --gM_flp-->                     0
                  0      0      0          0      -1     0        0         0          0       0     0            %Flp                      --gP_Flp-->                     0
                  0     -1      0          0       0     0        0         0          0       0     0            %Ara                      --gP_Ara-->                     0
                  0      0      1          0       0     0       -1         0          0       0     0            %DNA_1_BAD + Flp          --kR_Flp-->                     DNA_0_BAD + Flp
                  0     -1      0          0       0    -1        1         0          0       0     0            %DNA_1 + Ara              --k_BAD_on-->                   DNA_1_BAD
                  0      1      0          0       0     1       -1         0          0       0     0            %DNA_1_BAD                --k_BAD_off-->                  DNA_1 + Ara
                  0      0      0          0       0     0        0         1          1       0     0            %DNA_1_BAD                --kM_BAD_gfp_T7p-->             DNA_1_BAD + mRNA_gfp + mRNA_T7p
                  0      0      0          0       0     0        0         0          0       1     0            %mRNA_gfp                 --kP_GFP-->                     mRNA_gfp + GFP
                  0      0      0          0       0     0        0         0          0       0     1            %mRNA_T7p                 --kP_T7p-->                     mRNA_T7p + T7p
                  0      0      0          1       0     0        0         0          0       0     0            %DNA_1 + T7p              --kM_T7p_flp-->                 DNA_1 + T7p + mRNA_flp
                  0      0      0          1       0     0        0         0          0       0     0            %DNA_1_BAD + T7p          --kM_T7p_flp-->                 DNA_1_BAD + T7p + mRNA_flp 
                  1      0      0          0       0    -1        0         0          0       0     0            %DNA_1 + Flp              --kR_Flp-->                     DNA_0 + Flp
                 -1      0      0          0       0     1        0         0          0       0     0            %DNA_0 + Flp              --kR_Flp-->                     DNA_1 + Flp
                  0      0      0          0       0     0        0        -1          0       0     0            %mRNA_gfp                 --gM_gfp-->                     0
                  0      0      0          0       0     0        0         0         -1       0     0            %mRNA_T7p                 --gM_T7p-->                     0
                  0      0      0          0       0     0        0         0          0      -1     0            %GFP                      --gP_GFP-->                     0
                  0      0      0          0       0     0        0         0          0       0    -1];          %T7p                      --gP_T7p-->                     0

%% Run simulation
[t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
%[t,x] = firstReactionMethod(stoich_matrix, pfun, tspan, x0, p);


%% Plot time course
x_modified = x;

% REMOVE ARA COLUMN FROM VECTOR
% CAREFUL, INDEXES SHIFT WITH THIS OPERATION!
x_modified(:,2) = [];

% % CAUTION: ARA IS MISSING FROM THIS GRAPH
createfigure_1_0(t,x_modified, 'After second Ara pulse');

% CAUTION: only ARA being plotted
createfigure_ARA_degrade(t,x(:,2), 'Ara degradation after second pulse');

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
mRNA_T7p = x(9);
GFP = x(10);
T7p = x(11);

R_tot = p.k_BAD_on*DNA_0*Ara + ...
    p.k_BAD_off*DNA_0_BAD + ...
    p.kP_Flp*mRNA_flp + ...
    p.kR_Flp*DNA_0_BAD*Flp + ...
    p.gM_flp*mRNA_flp + ...
    p.gP_Flp*Flp + ...
    p.gP_Ara + ...
    p.kR_Flp*DNA_1_BAD*Flp + ...
    p.k_BAD_on*DNA_1*Ara + ...
    p.k_BAD_off*DNA_1_BAD + ...
    p.kM_BAD_gfp_T7p*DNA_1_BAD + ...
    p.kP_GFP*mRNA_gfp + ...
    p.kP_T7p*mRNA_T7p + ...
    p.kM_T7p_flp*DNA_1*T7p + ...  % is this even correct?
    p.kM_T7p_flp*DNA_1_BAD*T7p + ...
    p.kR_Flp*DNA_1*Flp + ...
    p.kR_Flp*DNA_0*Flp + ...
    p.gM_gfp*mRNA_gfp + ...
    p.gM_T7p*mRNA_T7p + ...
    p.gP_GFP*GFP + ...
    p.gP_T7p*T7p;

a = [p.k_BAD_on*DNA_0*Ara/R_tot;            %activation promoter
     p.k_BAD_off*DNA_0_BAD/R_tot;       %deactivation promoter
     p.kP_Flp*mRNA_flp/R_tot/R_tot;   %translation
     p.kR_Flp*DNA_0_BAD*Flp/R_tot;   %recombination
     p.gM_flp*mRNA_flp/R_tot;   %protein decay
     p.gP_Flp*Flp/R_tot; %protein decay 
     p.gP_Ara/R_tot; %protein decay !!! CONSTANT LIKE IN PAPER −cAra !!!
     %p.gP_Ara*Ara;   %protein decay !!! exponential −dAra ⋅[ara] !!!
     p.kR_Flp*DNA_1_BAD*Flp/R_tot;   %recombination
     p.k_BAD_on*DNA_1*Ara/R_tot; % activation promoter
     p.k_BAD_off*DNA_1_BAD/R_tot; % deactivation promoter
     p.kM_BAD_gfp_T7p*DNA_1_BAD/R_tot; 
     p.kP_GFP*mRNA_gfp/R_tot; 
     p.kP_T7p*mRNA_T7p/R_tot; 
     p.kM_T7p_flp*DNA_1*T7p/R_tot;
     p.kM_T7p_flp*DNA_1_BAD*T7p/R_tot; 
     p.kR_Flp*DNA_1*Flp/R_tot; 
     p.kR_Flp*DNA_0*Flp/R_tot; 
     p.gM_gfp*mRNA_gfp/R_tot; 
     p.gM_T7p*mRNA_T7p/R_tot; 
     p.gP_GFP*GFP/R_tot; 
     p.gP_T7p*T7p/R_tot;
     ]; 
end
