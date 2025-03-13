function [mag_out_rd, peak_mag, peak_idx]= RSP(fd, phase, nfft, Nsamp_nu, fft_prn, in)
sgl_dp_zp=zeros(nfft,1);
out_rd=zeros(nfft,length(fd));
% Multiplication of input signals with Doppler phases
sgl_dp=in.*phase;
% Matched Filtering one-by-one  for each Doppler frequency
for f=1:length(fd)
    sgl_dp_zp(1:Nsamp_nu,f)=sgl_dp(:,f);
    fft_sgl_dp(:,f)=fft(sgl_dp_zp(:,f),nfft);%FFT of received signal
   out_rd(:,f)=ifft(fft_sgl_dp(:,f).*fft_prn,nfft);  
end
%Peak Search for range and Doppler
mag_out_rd=abs(out_rd);
[peak_mag,peak_idx] = max(max(mag_out_rd));
end