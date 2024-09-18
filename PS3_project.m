A = 'Enter the desired method to solve the problem (NR - GS)';
disp(A)
method = input('Method: ', 's');  % 's' flag ensures that the input is treated as a string

if strcmp(method, 'NR')
         %% Ybus FORMULATION
Y(1,1) = 3.81563 -19.07814i +5.1696 -25.8478i +0.05125i +0.03875i;
Y(1,2) = -3.81563 +19.07814i;
Y(1,3)= -5.1696 +25.8478i;
Y(1,4) = 0;
Y(2,1) = -3.81563 +19.07814i;
Y(2,2) = 3.81563 -19.07814i + 5.1696 -25.8478i + 0.03875i +0.05125i;
Y(2,3) = 0;
Y(2,4) = -5.1696 +25.8478i;
Y(3,1) = -5.1696 +25.8478i;
Y(3,2) = 0; 
Y(3,3)= 5.1696 -25.8478i + 3.023705 -15.18528i + 0.06375i +0.03875i;
Y(3,4) =  -3.023705 +15.18528i;
Y(4,1) = 0;
Y(4,2) = -5.1696 +25.8478i;
Y(4,3) = -3.023705 +15.18528i;
Y(4,4) = 5.1696 -25.8478i +3.023705 -15.18528i + 0.06375i +0.03875i;
%% GIVENS
%slack bus:
V(1) =1+0i;
Slack_Bus_Load = 0.5 +0.3099i;
%PQ Buses:
Psch(2) = -1.7;
Qsch(2) = -1.0535;
V(2) = 1+0i; %Initial condition
Psch(3) = -2;
Qsch(3) = -1.2394;
V(3) = 1+0i; %Initial condition
%PV bus:
Psch(4) = 2.38;
V(4) = 1.02 + 0i;
%% ITERATIONS

for j = 0:1:5
    %Jacobian matrix calculation:
    
      jacobian(1,1) = abs(V(2)*Y(2,1)*V(1))*sin(angle(Y(2,1))+angle(V(1))-angle(V(2)))+abs(Y(2,4)*V(4)*V(2))*sin(angle(Y(2,4))+angle(V(4))-angle(V(2)));
      jacobian(1,2) = 0;
      jacobian(1,3) = -abs(V(2)*Y(2,4)*V(4))*sin(angle(Y(2,4))+angle(V(4))-angle(V(2)));
      jacobian(1,4) = 2*(abs(V(2))^2)*real(Y(2,2)) + abs(V(2)*Y(2,1)*V(1))*cos(angle(Y(2,1))+angle(V(1))-angle(V(2)))+abs(Y(2,4)*V(4)*V(2))*cos(angle(Y(2,4))+angle(V(4))-angle(V(2)));
      jacobian(1,5) =0;
      jacobian(2,1) =0;
      jacobian(2,2) = abs(V(3)*Y(3,1)*V(1))*sin(angle(Y(3,1))+angle(V(1))-angle(V(3)))+abs(Y(3,4)*V(4)*V(3))*sin(angle(Y(3,4))+angle(V(4))-angle(V(3)));
      jacobian(2,3) = -abs(V(3)*Y(3,4)*V(4))*sin(angle(Y(3,4))+angle(V(4))-angle(V(3)));
      jacobian(2,4) = 0;
      jacobian(2,5) = 2*(abs(V(3))^2)*real(Y(3,3))+ abs(V(3)*Y(3,1)*V(1))*cos(angle(Y(3,1))+angle(V(1))-angle(V(3)))+abs(Y(3,4)*V(4)*V(3))*cos(angle(Y(3,4))+angle(V(4))-angle(V(3)));
      jacobian(3,1) = -abs(V(2)*Y(4,2)*V(4))*sin(angle(Y(4,2))+angle(V(2))-angle(V(4)));
      jacobian(3,2) = -abs(V(3)*Y(4,3)*V(4))*sin(angle(Y(4,3))+angle(V(3))-angle(V(4)));
      jacobian(3,3) = abs(V(2)*Y(4,2)*V(4))*sin(angle(Y(4,2))+angle(V(2))-angle(V(4)))+abs(Y(4,3)*V(3)*V(4))*sin(angle(Y(4,3))+angle(V(3))-angle(V(4)));
      jacobian(3,4) = abs(V(4)*Y(4,2)*V(2))*cos(angle(Y(4,2))+angle(V(2))-angle(V(4)));
      jacobian(3,5) = abs(V(4)*Y(4,3)*V(3))*cos(angle(Y(4,3))+angle(V(3))-angle(V(4)));
      jacobian(4,1) = abs(V(2)*Y(2,1)*V(1))*cos(angle(Y(2,1))+angle(V(1))-angle(V(2)))+abs(Y(2,4)*V(4)*V(2))*cos(angle(Y(2,4))+angle(V(4))-angle(V(2)));
      jacobian(4,2) = 0;
      jacobian(4,3) =  -abs(V(2)*Y(2,4)*V(4))*cos(angle(Y(2,4))+angle(V(4))-angle(V(2)));
      jacobian(4,4) = -2*(abs(V(2))^2)*imag(Y(2,2)) -abs(V(2)*Y(2,1)*V(1))*sin(angle(Y(2,1))+angle(V(1))-angle(V(2)))-abs(Y(2,4)*V(4)*V(2))*sin(angle(Y(2,4))+angle(V(4))-angle(V(2)));
      jacobian(4,5) = 0;
      jacobian(5,1) = 0; 
      jacobian(5,2) = abs(V(3)*Y(3,1)*V(1))*cos(angle(Y(3,1))+angle(V(1))-angle(V(3)))+abs(Y(3,4)*V(4)*V(3))*cos(angle(Y(3,4))+angle(V(4))-angle(V(3)));
      jacobian(5,3) = -abs(V(3)*Y(3,4)*V(4))*cos(angle(Y(3,4))+angle(V(4))-angle(V(3)));
      jacobian(5,4) = 0;
      jacobian(5,5) = -2*(abs(V(3))^2)*imag(Y(3,3)) -abs(V(3)*Y(3,1)*V(1))*sin(angle(Y(3,1))+angle(V(1))-angle(V(3)))-abs(Y(3,4)*V(4)*V(3))*sin(angle(Y(3,4))+angle(V(4))-angle(V(3)));

    %Inverse jacobian matrix:
    Inv_Jacobian = inv(jacobian);
    
    %Calculated power(P):
    for cntr=2:1:4
    Psum = 0;
    for n=1:1:4
        Psum = Psum + abs(Y(cntr,n)*V(n))*cos(angle(Y(cntr,n))+angle(V(n))-angle(V(cntr)));
    end
    Pcalc(cntr) = abs(V(cntr))*Psum;
    
    end
    %Calculated reactive power(Q):
     for cntr=2:1:3
    Qsum = 0;
    for n=1:1:4
        Qsum = Qsum + abs(Y(cntr,n)*V(n))*sin(angle(Y(cntr,n))+angle(V(n))-angle(V(cntr)));
    end
    Qcalc(cntr) = -abs(V(cntr))*Qsum;
    
     end
    %Delta matrix calculation:
    delta_matrix(1,1) = Psch(2)-Pcalc(2);
    delta_matrix(2,1) = Psch(3)-Pcalc(3);
    delta_matrix(3,1) = Psch(4)-Pcalc(4);
    delta_matrix(4,1) = Qsch(2)-Qcalc(2);
    delta_matrix(5,1) = Qsch(3)-Qcalc(3);
    
    %Delta variables matrix:
    delta_var = Inv_Jacobian*delta_matrix;
    
    %New values:
    V2angle = angle(V(2))+ delta_var(1,1);
    V3angle = angle(V(3))+ delta_var(2,1);
    V4angle = angle(V(4))+ delta_var(3,1);
    V2mag = abs(V(2))*(1+ delta_var(4,1));
    V3mag = abs(V(3))*(1+ delta_var(5,1));
  
   
   V(2) = V2mag*(cos(V2angle)+i*sin(V2angle));
   V(3) = V3mag*(cos(V3angle)+i*sin(V3angle));
   V(4) = 1.02*(cos(V4angle)+i*sin(V4angle));
    
end

%% SLACK BUS POWER CALCULATION

CurrentSum =0;
for n=1:1:4
    CurrentSum = CurrentSum + V(n)*Y(1,n);
    
end
Slack_Bus_Injected_Power = V(1)*conj(CurrentSum);
Slack_Bus_Generated_Power = Slack_Bus_Injected_Power + Slack_Bus_Load ;

%% POWER FLOWS AND POWER LOSSES CALCULATIONS

%Line 1-2:
S12 = V(1)*conj(-(V(1)-V(2))*Y(1,2)+V(1)*0.05125i);
S21 = V(2)*conj(-(V(2)-V(1))*Y(1,2)+V(2)*0.05125i);
Line1_2Losses = S12+S21;
%Line 1-3:
S13 = V(1)*conj(-(V(1)-V(3))*Y(1,3)+V(1)*0.03875i);
S31 = V(3)*conj(-(V(3)-V(1))*Y(1,3)+V(3)*0.03875i);
Line1_3Losses = S13+S31;
%Line 3-4:
S34 = V(3)*conj(-(V(3)-V(4))*Y(3,4)+V(3)*0.06375i);
S43 = V(4)*conj(-(V(4)-V(3))*Y(3,4)+V(4)*0.06375i);
Line4_3Losses = S34+S43;
%Line 2-4:
S24 = V(2)*conj(-(V(2)-V(4))*Y(2,4)+V(2)*0.03875i);
S42 = V(4)*conj(-(V(4)-V(2))*Y(2,4)+V(4)*0.03875i);
Line4_2Losses = S24 +S42;

%% RESULTS DISPLAY

%Bus voltages display:
BUS_NUMBER = ["Bus 1";"Bus 2"; "Bus 3"; "Bus 4"];
VOLTAGE_Kv =[V(1)*230;V(2)*230;V(3)*230;V(4)*230];
table(BUS_NUMBER,VOLTAGE_Kv)

%Slack bus power display:
SLACK_BUS_POWER =["Injected power (MVA)";"Generated power (MVA)"];
S=[Slack_Bus_Injected_Power*100; Slack_Bus_Generated_Power*100];
table(SLACK_BUS_POWER,S)

%Line flows and lines losses display:
LINE_NUMBER =["Line 1-2";"Line 1-3";"Line 4-3";"Line 4-2";];
LINE_FLOW_MVA = [S12*100;S13*100;S43*100;S42*100];
DIRECTION_OF_FLOW = ["Bus 1 to Bus 2";"Bus 1 to Bus 3";"Bus 4 to Bus 3";"Bus 4 to Bus 2";];
LINE_LOSSES_MVA = [Line1_2Losses*100;Line1_3Losses*100;Line4_3Losses*100;Line4_2Losses*100];
table(LINE_NUMBER,LINE_FLOW_MVA,DIRECTION_OF_FLOW,LINE_LOSSES_MVA)

elseif strcmp(method, 'GS')
         %% Power flow analysis using Gauss method
Sbase = 100; Vbase = 230;
%% Lines data: Shunt Y is nominated by Yshunti_j
% Line 1-2:
R1_2 = 0.01008; X1_2 = 0.05040i; Yshunt1_2 = 0.05125i;
Y1_2 = 1/(R1_2 + X1_2);
% Line 1-3:
R1_3 = 0.00744; X1_3 = 0.03720i; Yshunt1_3 = 0.03975i;
Y1_3 = 1/(R1_3 + X1_3);
% Line 2-4:
R2_4 = 0.00744; X2_4 = 0.03720i; Yshunt2_4 = 0.03875i;
Y2_4 = 1/(R2_4 + X2_4);
% Line 3_4:
R3_4 = 0.01272; X3_4 = 0.06360i; Yshunt3_4 = 0.06375i;
Y3_4 = 1/(R3_4 + X3_4);
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Buses data:
% Bus 1: Slack bus
Pload_1 = 50; Qload_1 = 30.99; V_1 = 1;
Ploadpu_1 = Pload_1/Sbase; Qloadpu_1 = Qload_1/Sbase;
% Bus 2: Load bus (inductive)
Pload_2 = 170; Qload_2 = 105.35; V_2 = 1;
Ploadpu_2 = -Pload_2/Sbase; Qloadpu_2 = -Qload_2/Sbase;
% Bus 3: Load bus (inductive)
Pload_3 = 200; Qload_3 = 123.94; V_3 = 1;
Ploadpu_3 = -Pload_3/Sbase; Qloadpu_3 = -Qload_3/Sbase;
% Bus 4: Voltage controlled
Pgen_4 = 318; Pload_4 = 80; Qload_4 = 49.58; V_4 = 1.02;
Ploadpu_4 = (Pgen_4 - Pload_4)/Sbase; Q_4 = 0;
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Admittance bus:
Y11 = Y1_2 + Y1_3 + Yshunt1_2 + Yshunt1_3; Y12 = -Y1_2; Y13 = -Y1_3; Y14 = 0;
Y21 = -Y1_2; Y22 = Y1_2 + Y2_4 + Yshunt1_2 + Yshunt2_4; Y23 = 0; Y24 = -Y2_4;
Y31 = -Y1_3; Y32 = 0; Y33 = Y1_3 + Y3_4 + Yshunt1_3 + Yshunt3_4; Y34 = -Y3_4;
Y41 = 0; Y42 = -Y2_4; Y43 = -Y3_4; Y44 = Y2_4 + Y3_4 + Yshunt2_4 + Yshunt3_4;
Y = [ Y11 Y12 Y13 Y14 ; Y21 Y22 Y23 Y24 ; Y31 Y32 Y33 Y34 ; Y41 Y42 Y43 Y44];
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Define Vitertion array which contains the current iteration of each bus
%% voltage, which starts with the initial conditions as the following:
Viteration = [ V_1 ; V_2 ; V_3 ; V_4 ];
%% Define number of iterations: (15 iterations was selected to be the
%% required number of iterations till steady state results by trial & error)
Niteration = 15;
%% Bus voltage calculations based on the selected number of iterations:
for i = 1:Niteration
    for n = 2:4
         if n == 2
             Viteration(2) = (1/Y(2,2))*(((Ploadpu_2 - Qloadpu_2*1i)/(conj(Viteration(2)))) - Y(2,1)*Viteration(1) - Y(2,3)*Viteration(3) - Y(2,4)*Viteration(4));
         elseif n == 3
             Viteration(3) = (1/Y(3,3))*(((Ploadpu_3 - Qloadpu_3*1i)/(conj(Viteration(3)))) - Y(3,1)*Viteration(1) - Y(3,2)*Viteration(2) - Y(3,4)*Viteration(4));
         elseif n == 4
             Q_4 = -imag((conj(Viteration(4))*(Y(4,1)*Viteration(1) + Y(4,2)*Viteration(2) + Y(4,3)*Viteration(3) + Y(4,4)*Viteration(4))));
             Viteration(4) = (1/Y(4,4))*(((Ploadpu_4 - Q_4*1i)/(conj(Viteration(4)))) - Y(4,1)*Viteration(1) - Y(4,2)*Viteration(2) - Y(4,3)*Viteration(3));
             Viteration(4) = V_4*(cos(angle(Viteration(4))) + sin(angle(Viteration(4)))*1i);
        end
    end
end 
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Slack bus power calculation:
Pslack_inj = real(V_1*conj((Y(1,1)*V_1 + Y(1,2)*Viteration(2) + Y(1,3)*Viteration(3) + Y(1,4)*Viteration(4))));
Qslack_inj = imag(V_1*conj((Y(1,1)*V_1 + Y(1,2)*Viteration(2) + Y(1,3)*Viteration(3) + Y(1,4)*Viteration(4))));
Pslack_gen = Pslack_inj + Ploadpu_1;
Qslack_gen = Qslack_inj + Qloadpu_1;
Sslack_inj = Pslack_inj + Qslack_inj*1i;
Sslack_gen = Pslack_gen + Qslack_gen*1i;
%% Line flow & power losses:
% Line 1-2:
S1_2 = V_1*(conj((V_1 - Viteration(2))*Y1_2 + V_1*Yshunt1_2));
S2_1 = Viteration(2)*(conj((Viteration(2) - V_1)*Y1_2 + Viteration(2)*Yshunt1_2));
Slosses1_2 = S1_2 + S2_1;
% Line 1-3:
S1_3 = V_1*(conj((V_1 - Viteration(3))*Y1_3 + V_1*Yshunt1_3));
S3_1 = Viteration(3)*(conj((Viteration(3) - V_1)*Y1_3 + Viteration(3)*Yshunt1_3));
Slosses1_3 = S1_3 + S3_1;
% Line 2-4:
S2_4 = Viteration(2)*(conj((Viteration(2) - Viteration(4))*Y2_4 + Viteration(2)*Yshunt2_4));
S4_2 = Viteration(4)*(conj((Viteration(4) - Viteration(2))*Y2_4 + Viteration(4)*Yshunt2_4));
Slosses2_4 = S2_4 + S4_2;
% Line 3-4:
S3_4 = Viteration(3)*(conj((Viteration(3) - Viteration(4))*Y3_4 + Viteration(3)*Yshunt3_4));
S4_3 = Viteration(4)*(conj((Viteration(4) - Viteration(3))*Y3_4 + Viteration(4)*Yshunt3_4));
Slosses3_4 = S3_4 + S4_3;
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Results display:

%Bus voltages display:
Bus_number = ["Bus 1";"Bus 2"; "Bus 3"; "Bus 4"];
Bus_voltage_in_kV =Vbase*Viteration;
table(Bus_number,Bus_voltage_in_kV)

%Slack bus power display:
Slack_bus_data = ["Injected power in MVA";"Generated power in MVA"];
Slack_bus_power = Sbase*[Sslack_inj; Sslack_gen];
table(Slack_bus_data,Slack_bus_power)

%Line flows and lines losses display:
Line =["Line 1-2";"Line 1-3";"Line 2-4";"Line 3-4";];
Line_flow = Sbase*[S1_2; S1_3; S4_2; S4_3];
Direction = ["Bus 1 to Bus 2";"Bus 1 to Bus 3";"Bus 4 to Bus 3";"Bus 4 to Bus 2";];
Line_losses = Sbase*[Slosses1_2; Slosses1_3; Slosses2_4; Slosses3_4];
table(Line, Line_flow,Direction ,Line_losses)
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------

else
    disp('Invalid method. Please enter NR or GS.');
end
