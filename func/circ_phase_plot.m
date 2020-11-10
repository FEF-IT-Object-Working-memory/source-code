%% plot phase in circle
%% By Ehsan Rezayat

function hc = circ_phase_plot(theta)

hold on
polar([theta;theta],[zeros(size(theta));ones(size(theta))],'b')
                set(gca,'box','off')
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                text(1.2, 0, '0'); text(-.05, 1.2, '\pi/2');  text(-1.35, 0, '±\pi');  text(-.075, -1.2, '-\pi/2');
                axis square;
                zz = exp(1i*linspace(0, 2*pi, 101));                
                plot( real(zz), imag(zz), 'k', [-2 2], [0 0], 'k:', [0 0], [-2 2], 'k:');
                set(gca, 'XLim', [-1.1 1.1], 'YLim', [-1.1 1.1])
                set(gca, 'Visible', 'off')
end
  