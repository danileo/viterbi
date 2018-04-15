close all;
clear;

mex viterbi.c

% -- simulation parameters -- %
random_message=0; % 1 for random input message, 0 for all -1 bits (useful for repeteability)
Gamma_dB=-2:1:15; % SNR range in dB
Gamma=10.^(Gamma_dB/10); % SNR
k=1e4; % packet length
Nit=1e5; % iterations number
Th_nerr=1e3; % do not simulate if Th_nerr errors are exceeded

% -- code parameters -- %
global deg;
deg=2; % = nu
window_size=10; % window size (for example, window_size=10=5*deg if deg=2)
n=2; % = 1/R
p1=3; % \ puncture parameters such that
p2=4; % / the overall code rate is p1/p2
p=gcd(n*p1,p2);
q1=(n*p1)/p;
q2=p2/p;

prec=1+log10(1/(k*Nit)); % expected precision
disp(['Expected precision = 10^' num2str(prec)])

% initialize error and packet counter
nerr=zeros(size(Gamma));
nerr_punct=zeros(size(Gamma));
nerr_window=zeros(size(Gamma));
npack=zeros(size(Gamma));

tic % start time counter

% -- main cycle --
for it=1:Nit
    if random_message
        % randomly build binary message
        u=randi(2,1,k)-1;
        y=code(u);
        % map with map L
        s=2*y-1;
    else
        s=-ones(1,2*(k+deg));
    end
    
    % prepare noise samples (with unit variance)
    w=randn(size(s));
    
    % -- cycle on SNRs --
    for m=1:length(Gamma)
        
        % check number of errors
        if nerr(m)<Th_nerr
            
            % define noise variance
            sigw=sqrt(1/Gamma(m));
            % define the received signal
            r=s+sigw*w;
            
            % puncture code
            if p1>0 && p2>0
                r_punc=zeros(1,floor(length(r)*p2/(n*p1)));
                q=1;
                for p=1:length(r)
                    if mod(p-1,q1)<q2
                        r_punc(q)=r(p);
                        q=q+1;
                    end
                end
            end
            
            % demodulate
            if window_size>0
                [u_hat u_hat_window]=viterbi(r,[0 0],window_size);
            else
                u_hat=viterbi(r,[0 0],0);
            end
            if p1>0 && p2>0
                u_hat_punct=viterbi(r_punc,[p1 p2],0);
            end
            
            % count errors
            %nerr(m)=nerr(m)+sum(u~=u_hat);
            nerr(m)=nerr(m)+sum(u_hat);
            if p1>0 && p2>0
                %l=min(length(u),length(u_hat_punct));
                %nerr_punct(m)=nerr_punct(m)+sum(u(1:l)~=u_hat_punct(1:l));
                nerr_punct(m)=nerr_punct(m)+sum(u_hat_punct);
            end
            if window_size>0
                %nerr_window(m)=nerr_window(m)+sum(u~=u_hat);
                nerr_window(m)=nerr_window(m)+sum(u_hat_window);
            end
            npack(m)=npack(m)+1;
            
        end
        
    end
    
    if mod(it,100) == 0
        disp(['#' num2str(it) ', BER = ' num2str(nerr./(npack*k))]);
        save('uncoded');
    end
end

toc % read time counter

% calculate BER
Pbit=nerr./(npack*k);
if p1>0 && p2>0
    Pbit_punct=nerr_punct./(npack*k);
else
    Pbit_punct=zeros(1,length(Gamma_dB));
    p2=0;
    p1=1;
end
if window_size>0
    Pbit_window=nerr_window./(npack*k);
else
    Pbit_window=zeros(1,length(Gamma_dB));
end

% expected BER (uncoded)
Q=@(x) 0.5*erfc((x)/sqrt(2));
Pexp_uncoded=Q(sqrt(Gamma));

% expected BER (coded)
Pexp_coded=(3-6/k)*Q(sqrt(5*Gamma));

% show results
figure;
set(0,'defaultTextInterpreter','latex')
set(gca,'FontSize',14);

semilogy(Gamma_dB-10*log10(2),Pexp_uncoded,'--','Color',[.7 0 0])
hold on;
semilogy(Gamma_dB-10*log10(2/n),Pexp_coded,'k--')
semilogy(Gamma_dB-10*log10(2/n),Pbit,'k-')
semilogy(Gamma_dB-10*log10(2*p1/p2),Pbit_punct,'b-')
semilogy(Gamma_dB-10*log10(2/n),Pbit_window,'m*','MarkerSize',8)
semilogy(Gamma_dB-10*log10(2/n),ones(size(Gamma_dB))*10^prec,'r--')

axis([min(Gamma_dB-10*log10(2/n)) max(Gamma_dB-10*log10(2/n)) 1e-7 1e0])
legend('Uncoded (closed form)','Coded (closed form)', ...
    'Coded (simulated)','Coded, punctured version (simulated)', ...
    'Coded, windowed version (simulated)')
xlabel('$E_b/N_0$ [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca,'XMinorTick','on','YMinorTick','on','Ygrid','on','XGrid','on', ...
    'xcolor',[.5 .5 .5],'ycolor',[.5 .5 .5]);
Caxes = copyobj(gca,gcf);
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 38 15])
print -depsc plot.eps
