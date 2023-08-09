
c = 1:32; % initial concentration gradient (mol/(m^2*s))
D = 2;    % diffusivity (m^2/s)
dt = 0.1; % time step (s)
dx = 1;   % x step (m)
        
for k = 0:1000,
  % finite-difference diffusion using explicit method, central differences
  c = c + D*dt/(dx^2)*([c(2:end) c(end)] - 2*c + [c(1) c(1:end-1)]);

  % The MATLAB *plotting* code presented in the lesson does not work inside the Jupyter-notebook
  % Octave environment. Instead of using the "image" command, we need to "fill" individual boxes. 
  vertices = [0 0; 1 0; 1 1; 0 1];
  if mod(k,100) == 0, % plot a snapshot
    subplot(11,1,k/100+1); 
    for m = 1:length(c),
      fill(vertices(:,1)+m, vertices(:,2)-k/100,c(m)); hold on
    end
    caxis([1 32])  
    h = ylabel(sprintf('t = %gs      ',k*dt)); 
    set(gca,'ytick',[],'xticklabel',[],'ticklength',[0 0]); grid on
    set(gca,'xtick',1.5:1:100,'gridlinestyle','-','linewidth',4);
    set(h,'rotation',0,'horizontalalignment','right','verticalalignment','middle')
  end
end
xlabel('x location'); text(16,-14.25,'Diffusion example','horizontalalignment','center');

max(c) - min(c)
