function MI=MI_(phase_bias,power_bias)

  n_hist_bins = 18;
            phase_edges=linspace(-pi,pi,n_hist_bins+1);
            amp_by_phases=zeros(1,n_hist_bins);
            for i=1:n_hist_bins-1
                amp_by_phases(i) = nanmean(power_bias(phase_bias>phase_edges(i) & phase_bias<phase_edges(i+1)));
            end
            P_r=amp_by_phases/nansum(amp_by_phases);
            P_r=P_r(1:end-1);
            MI=(1*nansum(P_r.*log(P_r))+log(length(P_r)))/log(length(P_r));
end