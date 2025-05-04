clear all
clc

% Read the content
inputFileName = 'Experiment_Adaptive_DMC_P80.csv';
fileID = fopen(inputFileName, 'r');
if fileID == -1
    error('Cannot open input file');
end
fileContent = fread(fileID, '*char')';
fclose(fileID);
% Change ',' with '.'
modifiedContent = strrep(fileContent, ',', '.');
% Write the content
fileID = fopen(inputFileName, 'w');
if fileID == -1
    error('Cannot open output file');
end
fwrite(fileID, modifiedContent, 'char');
fclose(fileID);

Data = readtable(inputFileName,'NumHeaderLines',1);

[size_x, ~] = size(Data);

T = [];
T_IN = [];
T_SP = [];
OR_DMC = [];

for i=1:size_x
    if(strcmp(Data.Var1(i), 'Temperature_Out'))
        T = [T; Data.Var3(i)];
    elseif(strcmp(Data.Var1(i), 'Temperature_In'))
        T_IN = [T_IN; Data.Var3(i)];
    elseif(strcmp(Data.Var1(i), 'Temperature_Setpoint'))
        T_SP = [T_SP; Data.Var3(i)];
    elseif(strcmp(Data.Var1(i), 'OR'))
        OR_DMC = [OR_DMC; Data.Var3(i)];
    end
end

[size_xT, ~] = size(T);
[size_xT_IN, ~] = size(T_IN);
[size_xT_SP, ~] = size(T_SP);
[size_xOR, ~] = size(OR_DMC);

zeroLogicalArray = (OR_DMC == OR_DMC(1));
zerosNumber = sum(zeroLogicalArray);

if(zerosNumber == 0)
    zerosNumber = 1;
end

length = min([size_xT, size_xT_SP, size_xOR]);
T = T(zerosNumber:length,1);
T_SP = T_SP(zerosNumber:length,1);
OR_DMC = OR_DMC(zerosNumber:length,1);

DMC_offset = 20;
T = T(DMC_offset:end);
T_SP = T_SP(DMC_offset:end);
OR_DMC = OR_DMC(DMC_offset:end);

time = (0:1:length-zerosNumber-DMC_offset+1)';

T_OP = zeros(size(T,1), 1) + 50;

% Read the content
inputFileName = 'Experiment_PID_P80.csv';
fileID = fopen(inputFileName, 'r');
if fileID == -1
    error('Cannot open input file');
end
fileContent = fread(fileID, '*char')';
fclose(fileID);
% Change ',' with '.'
modifiedContent = strrep(fileContent, ',', '.');
% Write the content
fileID = fopen(inputFileName, 'w');
if fileID == -1
    error('Cannot open output file');
end
fwrite(fileID, modifiedContent, 'char');
fclose(fileID);

Data = readtable(inputFileName,'NumHeaderLines',1);

[size_x, ~] = size(Data);

T_PID = [];
T_IN_PID = [];
T_SP_PID = [];
Flow_PID = [];

for i=1:size_x
    if(strcmp(Data.Var1(i), 'rTemperatureOut'))
        T_PID = [T_PID; Data.Var3(i)];
    elseif(strcmp(Data.Var1(i), 'rTemperatureIn'))
        T_IN_PID = [T_IN_PID; Data.Var3(i)];
    elseif(strcmp(Data.Var1(i), 'rTemperatureSetpoint'))
        T_SP_PID = [T_SP_PID; Data.Var3(i)];
    elseif(strcmp(Data.Var1(i), 'rFlow'))
        Flow_PID = [Flow_PID; Data.Var3(i)];
    end
end

[size_xT, ~] = size(T_PID);
[size_xT_IN, ~] = size(T_IN_PID);
[size_xT_SP, ~] = size(T_SP_PID);
[size_xFlow, ~] = size(Flow_PID);

zeroLogicalArray = (Flow_PID == Flow_PID(1));
zerosNumber = sum(zeroLogicalArray);

if(zerosNumber == 0)
    zerosNumber = 1;
end

length = min([size_xT, size_xT_SP, size_xFlow]);
T_PID = T_PID(zerosNumber:length,1);
T_SP_PID = T_SP_PID(zerosNumber:length,1);
Flow_PID = Flow_PID(zerosNumber:length,1);

PID_offset = 1;

T_PID = T_PID(PID_offset:end);
T_SP_PID = T_SP_PID(PID_offset:end);
Flow_PID = Flow_PID(PID_offset:end);

time_PID = (0:1:length-zerosNumber-PID_offset+1)';

T_OP = zeros(size(T,1), 1) + 50;

f = figure;
f.Position = [100 100 1200 900];

subplot(2,1,1)
stairs(time, T_SP, '--', 'LineWidth', 1.1, 'Color', '#FF7F00');
hold on
plot(time_PID, T_PID, 'LineWidth', 1.5, 'Color', '#34831B');
plot(time, T, 'LineWidth', 1.5, 'Color', 'b');

% yline(50+0.05*5, 'r')
% yline(50-0.05*5, 'r')

OR_PID = Flow_PID*100/12;

hold off
% xlim([0 max(time)])
xlim([0 5600])
% ylim([44 51])
xlabel('t, s')
yl = sprintf("T, %cC", char(176));
ylabel(yl);
legend('Punkt pracy', 'Regulacja PID', 'Adaptacyjna regulacja DMC', 'Location', 'southeast')
grid on
subplot(2,1,2)
stairs(time_PID, OR_PID, 'LineWidth', 1.5, 'Color', '#34831B');
hold on
stairs(time, OR_DMC, 'LineWidth', 1.5, 'Color', 'b');
hold off
% xlim([0 max(time)])
xlim([0 5600])
xlabel("t, s");
ylabel("OR, %");
legend('Sterowanie - regulator PID', 'Sterowanie - adaptacyjny regulator DMC', 'Location', 'northeast')
grid on

%% Criteria
clc
% u
disp("max u DMC: " + num2str(max(P(:))))
disp("max u PID: " + num2str(max(Flow_PID(:))))

% du
diff_v = diff(P);
diff_v_PID = diff(Flow_PID);
max_change = max(abs(diff_v));
max_change_PID = max(abs(diff_v_PID));
disp("max du DMC: " + num2str(max_change))
disp("max du PID: " + num2str(max_change_PID))

% IAE
IAE = 0;
IAE_PID = 0;
length_IAE = min(numel(T),numel(T_PID));
for i=1:length_IAE
    IAE = IAE + i*abs(T(i) - T_OP(1));
    IAE_PID = IAE_PID + i*abs(T_PID(i) - T_OP(1));
end
disp("IAE DMC: " + num2str(IAE))
disp("IAE PID: " + num2str(IAE_PID))

% Max overshoot
ov = (max(T) - 50)/5 * 100;
ov_PID = (max(T_PID) - 50)/5 * 100;
disp("Max overshoot DMC: " + num2str(ov) + "%")
disp("Max overshoot PID: " + num2str(ov_PID) + "%")
