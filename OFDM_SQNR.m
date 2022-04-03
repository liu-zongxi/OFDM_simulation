%-----------------------OFDM中SNQR的计算-------------------%
%-----------------------author:lzx-------------------------%
%-----------------------date:2022年3月31日-----------------%

%%
clear, clf, clc;
Nfft = 64;
Nk = 64;                % 这三个参数不重复说明啦
Npsk = 6;
M = 2^Npsk;
MaxIter=1000;           % 迭代次数
Quantbits = [6:9];      % 量化比特数
V_cutoffs = [2:0.2:8];  % 限幅电平
gss=['ko-';'ks-';'k^-';'kd-'];
%%
for ii = 1:length(Quantbits)
    Quantbit = Quantbits(ii);
    Quantbit_fractional = Quantbit -1;
    for jj = 1:length(V_cutoffs)
        V_cutoff = V_cutoffs(jj);
        Tx = 0; Te = 0;
        for kk = 1:MaxIter
            X_mod = ModSymbolGenerator(Npsk, Nk);
            x = ifft(X_mod, Nfft);
            x = x/sqrt(1/(2*Nk))/V_cutoff;          % 实部或者虚部的方差是(1/(2*Nk))
            % x(:) = (real(x(:))>V_cutoff).*V_cutoff;
%             for kkk= 1:length(x)
%                 if real(x(kkk))>V_cutoff
%                     x(kkk) = V_cutoff + 1j*imag(x(kkk));
%                 end
%                 if imag(x(kkk))>V_cutoff
%                     x(kkk) = real(x(kkk)) + 1j*V_cutoff;
%                 end
%             end
            xq=fi(x,1,Quantbit,Quantbit_fractional)
            xq=double(xq); Px = sum(abs(x).^2); e = x-xq; Pe = sum(abs(e).^2);
            Tx = Tx + Px;  Te = Te + Pe;
        end
        SQNRdB(ii,jj) = 10*log10(Tx/Te);
    end
end
%%
[SQNRdBmax,imax] = max(SQNRdB'); % To find the maximum elements in each row of SQNRdB
for i=1:size(gss,1), plot(V_cutoffs,SQNRdB(i,:),gss(i,:)), hold on;  end
for i=1:size(gss,1)
   str(i,:)=[num2str(Quantbits(i)) 'bit quantization'];
   plot(V_cutoffs(imax(i)),SQNRdBmax(i),gss(i,1:2),'markerfacecolor','r')
end
xlabel('mus(clipping level normalized to \sigma)');  ylabel('SQNR[dB]');
legend(str(1,:),str(2,:),str(3,:),str(4,:));  grid on