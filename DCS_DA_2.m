%% Signal Generation

% First part of the Coding Involves generating the values of SAR given the
% Specific Heat of the Body and the Rate of Temperature Change

clc;
clear all;
syms x y;
% Generating SAR Values with the formula
C = (3.5*10^3); % Specific Heat Capacity of Human Body (3.5 KJ/Kg K)
rate_change = [0, 1.773 , 0.012 , 0.370 , 0.251 , 0.523 , 4.443 , 1.505 , 1.115 , 0.082 , 1.464 , 1.059 , 0.524 , 1.126 , 0.027 , 3.330 ]*10^-4; % Rate Change has been computed every hour and with a distance of 100mm 
Average_SAR_Values = C*rate_change;
V = max(Average_SAR_Values);
t = [0:1:15]; % nth Hour
% Plotting the Curve of Time vs SAR
stem(t , Average_SAR_Values,'LineWidth' , 3, 'Color', [0.905, 0.517, 0.517])
ylabel('Specific Absorption Rate (Watt/Kg)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time(Hour)', 'FontSize', 12, 'FontWeight', 'bold');
title('SAR Values Plot ', 'FontSize', 12, 'FontWeight', 'bold');
set(gca,'FontSize',13,'FontWeight','bold');


%% Performing Non Uniform Quantization using mu Law

clc;
clear all;
syms x y;
% Generating SAR Values with the formula
C = (3.5*10^3); % Specific Heat Capacity of Human Body (3.5 KJ/Kg K)
rate_change = [0, 1.773 , 0.012 , 0.370 , 0.251 , 0.523 , 4.443 , 1.505 , 1.115 , 0.082 , 1.464 , 1.059 , 0.524 , 1.126 , 0.027 , 3.330 ]*10^-4; % Rate Change has been computed every hour and with a distance of 100mm 
Average_SAR_Values = C*rate_change;
V = max(Average_SAR_Values);
mu = 255;
partition = 0:0.01:1.7;
codebook = 0:1.7/171:1.7;
[~,~,distortion] = quantiz(Average_SAR_Values,partition,codebook); % distortion gives out the bias between the quantised and the original value
compsig = compand(Average_SAR_Values,mu,V,'mu/compressor'); % Compressing using mu law
[~,quants] = quantiz(compsig,partition,codebook);
newsig = compand(quants,mu,max(quants),'mu/expander'); % Expanding using mu law
distortion2 = sum((newsig - Average_SAR_Values).^2)/length(Average_SAR_Values);
stem([Average_SAR_Values' compsig']);
set(gca,'FontSize',5,'FontWeight','bold');
title('Comparison between Original signal and Companded Signal using mu law');
xlabel('Interval');
ylabel('Amplitude');
legend('Original','Companded','location','nw');

%% Performing Non Uniform Quantization using A law

clc;
clear all;
syms x y;
% Generating SAR Values with the formula
C = (3.5*10^3); % Specific Heat Capacity of Human Body (3.5 KJ/Kg K)
rate_change = [0, 1.773 , 0.012 , 0.370 , 0.251 , 0.523 , 4.443 , 1.505 , 1.115 , 0.082 , 1.464 , 1.059 , 0.524 , 1.126 , 0.027 , 3.330 ]*10^-4; % Rate Change has been computed every hour and with a distance of 100mm 
Average_SAR_Values = C*rate_change;
V = max(Average_SAR_Values);
A = 255;
partition = 0:0.01:1.7;
codebook = 0:171;
[~,~,distortion] = quantiz(Average_SAR_Values,partition,codebook);
compsig = compand(Average_SAR_Values,A,V,'A/compressor');
[~,quants] = quantiz(compsig,partition,codebook);
newsig = compand(quants,A,max(quants),'A/expander');
distortion2 = sum((newsig - Average_SAR_Values).^2)/length(Average_SAR_Values);
stem([Average_SAR_Values' compsig']);
title('Comparison between Original signal and Companded Signal using A Law');
xlabel('Interval');
ylabel('Apmlitude');
legend('Original','Companded','location','nw');

%% Performing Uniform Quantization

clc;
clear all;
syms x y;
% Generating SAR Values with the formula
C = (3.5*10^3); % Specific Heat Capacity of Human Body (3.5 KJ/Kg K)
rate_change = [0, 1.773 , 0.012 , 0.370 , 0.251 , 0.523 , 4.443 , 1.505 , 1.115 , 0.082 , 1.464 , 1.059 , 0.524 , 1.126 , 0.027 , 3.330 ]*10^-4; % Rate Change has been computed every hour and with a distance of 100mm 
Average_SAR_Values = C*rate_change;
V = max(Average_SAR_Values);
t = [0:1:15];
partition = 0:0.01:1.7;
codebook = 0:1.7/171:1.7;
[index,quants] = quantiz(Average_SAR_Values,partition,codebook);
stem([quants' , Average_SAR_Values'])
title('Quantization of SAR Signal ')
xlabel('Time')
ylabel('Amplitude')
legend('Original sampled SAR values','Quantized SAR Values');


% If we consider Error then uniform Quantization performs the best and if
% we consider complexity into consideration then we can choose non uniform
% quantizer but the major problem that is needed to be tackled in the case
% of non uniform quantizer is the error or bias introduced in the quantized
% signal compared to that of the original signal.


%% Modulation using PCM and Demodulation

clc;
clear all;
% Generating SAR Values with the formula
C = (3.5*10^3); % Specific Heat Capacity of Human Body (3.5 KJ/Kg K)
rate_change = [0, 1.773 , 0.012 , 0.370 , 0.251 , 0.523 , 4.443 , 1.505 , 1.115 , 0.082 , 1.464 , 1.059 , 0.524 , 1.126 , 0.027 , 3.330 ]*10^-4; % Rate Change has been computed every hour and with a distance of 100mm 
Average_SAR_Values = C*rate_change;
xmin = 0;
xmax = 1.7;
n = 5; % 5 bit Encoding or 5 bits per sample
L = 2^n; % No of quantization Levels
del = (xmax - xmin)/2;

V = max(Average_SAR_Values);
t = [0:1:15];
partition = 0:1.7/L:1.7;
codebook = 0:1.7/L:1.7+1.7/L;
[index,quants] = quantiz(Average_SAR_Values,partition,codebook);
l1=length(index);     			  % to convert 1 to n as 0 to n-1 indices
for i=1:l1
    if (index(i)~=0)
        index(i)=index(i)-1;
    end
end
l2=length(quants);
for i=1:l2 			%  to convert the end representation levels within the range.
    if(quants(i)==xmin-(del/2))
        quants(i)=xmin+(del/2);
    end
    if(quants(i)==xmax+(del/2))
        quants(i)=xmax-(del/2);
    end
end
code = de2bi(index,'left-msb') 	% Decimal to Binary Conversion
k=1;
for i=1:l1                     % to convert column vector to row vector
    for j=1:n
        coded(k)=code(i,j);
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
% Demodulation

code1=reshape(coded,n,(length(coded)/n));
index1=bi2de(code1,'left-msb');
re_signal=del*index+xmin+(del/2);
figure(5);
hold on;
stem([(re_signal*0.064)' Average_SAR_Values'])
title('Demodulated Signal');
xlabel('time');
ylabel('Amplitude');
hold off;

%% Delta Modulation

clc;
clear all;
close all;
rate_change = [0, 1.773 , 0.012 , 0.370 , 0.251 , 0.523 , 4.443 , 1.505 , 1.115 , 0.082 , 1.464 , 1.059 , 0.524 , 1.126 , 0.027 , 3.330 ]*10^-4; % Rate Change has been computed every hour and with a distance of 100mm 
C = (3.5*10^3);
Average_SAR_Values = C*rate_change;% input signal
l=length(Average_SAR_Values);
xmin = 0;
xmax = 1.7;
figure(1);
plot(Average_SAR_Values,'g');
title("Input Signal");
delta=xmax*0.01; % assigning delta value for step changes in staircase
xn=0; % initial assumption of staircase 
for i=1:l
if Average_SAR_Values(i)>xn(i) % conditional check between i/p and staircase
d(i)=1; xn(i+1)=xn(i)+delta; % transmit bit 1 and Add +delta with staircase
else
d(i)=0; xn(i+1)=xn(i)-delta; % transmit bit 0 and Add +delta with staircase
end
end
figure;
stairs(xn) % Approximation of x(t) in staircase
title("Delta Modulation");
y=[d] % transmitted bit sequence
figure;
stairs(y)
title("Bit Sequence Transmitted");
l2=length(y); 
xn1=0  % reconstructed staircase (intially assign zero)

% demodulation

for i=1:l2   
if y(i)==1
xn1(i+1)=xn1(i)+delta;
else
xn1(i+1)=xn1(i)-delta;
end
end
figure;
stairs(xn1)
title("Demodulation")
[num,den]=butter(2,0.1,'low');%demodulation using a low Pass Filter
Finale=filter(num,den,xn1); %filtering the signal
figure(5);
plot(Finale,'b');
title('Demodulation using Lowpass Filter');
