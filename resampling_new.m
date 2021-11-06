function sig=resampling_new(sig_or,f1,f2,num_chunks,tail_length,visualize)
% function sig=resampling(sig_or,f1,f2,num_chunks,tail_length,visualize)
%
% Resampling function based on interpolating data with an spline-like
% function and reconstructing signal at the new sampling frequency.
%
% sig_or: original signal
%
% f1: original sampling frequency
%
% f2: new sampling frequency
%
% num_chunks: data can be divided into portions (max 4 to avoid disalignments) 
%    to speed up computation and reduce memory requirements
%
% tail_length: the last n samples must be treated differently. Insert 0 to
%   use default value (20000) or increase it if a divergent signal appears
%
% visualize: insert 1 to plot the original and the resampled signals,
%   otherwise 0

if num_chunks>4
    error('Number of chunks must be lower than 4');
else

    sig=sig_or;
    t_or=1/f1:1/f1:length(sig_or)/f1;
    if tail_length==0
        tail_length=20000;
    end
    max_length=round((length(sig)-3*tail_length)/num_chunks);
    steps=floor(length(sig)/max_length);
    if mod(length(sig),max_length)>2*tail_length
        steps=steps+1;
    end
    index=1;
    sig_=[];
    for i=1:steps-1
        sig_{i}=sig(index:index+max_length);
        index=index+max_length-tail_length;
    end
    sig_{steps}=sig(index:length(sig));
    sig_end=sig(end-tail_length:end);
    for i=1:3
        sig_end=[sig_end,sig_end];
    end
    for i=1:length(sig_)
        sig_app=sig_{i};
        t=1/f1:1/f1:length(sig_app)/f1;
        pp=pchip(t,sig_app);
        t2=1/f2:1/f2:max(t);
        clear t;
        sig_{i}=ppval(pp,t2);
        clear t2;
    end
    t=1/f1:1/f1:length(sig_end)/f1;
    t2=1/f2:1/f2:max(t);
    pp=pchip(t,sig_end);
    sig_end=ppval(pp,t2);
    index=1;
    sig=[];
    for i=1:steps-1
        app=sig_{i};
        sig(index:index+round(max_length*f2/f1)-1)=app(1:round(max_length*f2/f1));
        index=index+round((max_length-tail_length)*f2/f1);
    end
    sig(index:index+length(sig_{steps})-1)=sig_{steps};
    sig(end-round(tail_length*f2/f1):end)=[];
    sig(end+1:end+round(tail_length*f2/f1))=sig_end(1:round(tail_length*f2/f1));
    
end

if visualize==1
    figure;
    plot(t_or,sig_or)
    hold on
    t_new=1/f2:1/f2:length(sig)/f2;
    plot(t_new,sig,'r');
    legend('Original signal','Resampled signal');
end
