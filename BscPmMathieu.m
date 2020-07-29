classdef BscPmMathieu < ott.BscPmParaxial
% Generate a focussed Mathieu beam using paraxial point matching
% Inherits from :cls:`ott.BscPmParaxial`.
%
% This code is described/used in:
%   CITATION
% If you find this code useful, please consider citing the paper.
%
% Properties:
%   - order         -- (numeric) Mathieu beam order
%   - parity        -- (enum) Mathieu beam parity (odd, even, or helical+/-)
%   - ellipticity   -- (numeric) Ellipticity of beam
%   - scale         -- (numeric) Scaling relative to lens aperture
%   - ifocal        -- (numeric) Interfocal spacing
%
% Properties (inherited):
%   - a             -- (complex) Beam shape coefficients a vector
%   - b             -- (complex) Beam shape coefficients b vector
%   - type          -- (enum) Beam type (incoming, outgoing or scattered)
%   - polarisation  -- (2xcomplex) Polarisation Jones vector
%   - NA            -- (numeric) Numerical aperture of focussing lens
%
% Static methods:
%
% See also BscPmMathieu, ott.BscPmParaxial, ott.BscPmGauss

% This file is part of the optical tweezers toolbox.
% Based on code by M. A Martinez-Ruiz INAOE, Puebla, Mex.,
% Isaac C.D. Lenton, The University of Queensland & Ulises-Ruiz INAOE
% Puebla, Mexico!!!
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    order         % (numeric) Mathieu beam order
    parity        % (enum) Mathieu beam parity (odd, even, or helical)
    ellipticity   % (numeric) Ellipticity of beam
    scale         % (numeric) Scaling relative to lens aperture
    ifocal        % (numeric) Inter-focal spacing
  end

  methods (Static)
    
    function [u, v] = cart2ellip(f, x, y)
      % Convert from Cartesian to elliptical coordinates
      %
      % Usage:
      %   [u, v] = cart2ellip(f, x, y)
      %
      % Parameters:
      %   f -- (numeric) separation between ellipse focal points
      %   x,y -- (numeric) Cartesian coordinates [-Inf, Inf]
      
      assert(all(size(x) == size(y)) || all(size(x) == 1 | size(y) == 1), ...
        'size of x and y must match or they must be transpose vectors');
      
      rho = sqrt(x.^2 + y.^2);
      
      u =acosh(sqrt((rho.^2 + f^2 + sqrt((rho.^2 + f^2).^2 - 4*(f^2)*(x.^2)))/(2*(f^2))));
      u=abs(u);
      
      v = zeros(size(rho));
      if any(size(x) ~= size(rho))
        x = repmat(x, numel(rho)./numel(x), 1);
        x = reshape(x, size(rho));
      end
      
      % second and third quadrant
      idx = x <= 0;
      v(idx) = acos(-sqrt((rho(idx).^2 + f^2 - sqrt((rho(idx).^2 + f^2).^2 - 4*(f^2)*(x(idx).^2)))/(2*(f^2))));
      
      % first and fourth quadrant
      idx = x > 0;
      v(idx) = acos(sqrt((rho(idx).^2 + f^2 - sqrt((rho(idx).^2 + f^2).^2 - 4*(f^2)*(x(idx).^2)))/(2*(f^2))));
      
      v = v .* sign(y);
    end
    
    function [x,y] = ellip2cart(f, u, v)
      % Convert from elliptical coordinates to Cartesian coordiantes
      %
      % Usage:
      %   [x, y] = ellip2cart(f, u, v)
      %
      % Parameters:
      %   f -- (numeric) separation between ellipse focal points
      %   u -- (numeric) elliptic coordinate [0, Inf]
      %   v -- (numeric) elliptic coordinate [-pi, pi] or equiv.

      x = f*cosh(u).*cos(v);
      y = f*sinh(u).*sin(v);

    end

    function eScalar = calculate_scalar_field(...
        parity, order, nscale, ellipticity, ifocal, varargin)
      % Calculate the scalar field for the Mathieu beam
      %
      % Usage:
      %   eScalar = calculate_scalar_field(parity, order,
      %   nscale, ellipticity, ifocal)
      %
      % Parameters:
      %   - order -- (numeric) Mathieu beam order
      %   - parity -- (enum) Mathieu beam parity (odd, even or helical+/-)
      %   - scale -- (numeric) Scaling factor for pattern.
      %   - ifocal -- (numeric) Interfocal separation [0, 1].
      %
      %   - ellipticity -- (numeric) Ellipticity of beam (q).
      %     This is related to the interfocal separation and the transverse
      %     spatial frequency by `q = f^2 k_T^2 / 4` where f is the semi-focal
      %     distance and k_T is the transverse spatial frequency.
      %
      %   - Nres -- Resolution for paraxial field calculation.
      %     Default: 256.
      
      p = inputParser;
      p.addParameter('Nres', 256);
      p.addParameter('projection', 'linear');
      p.parse(varargin{:});
      
      k = 2*pi;
      kT = sqrt(ellipticity).*2./ifocal;
      kL = k.^2 - kT.^2;
      
      assert(ellipticity >= 0, 'Ellipticity must be positive');
%       assert(ifocal <= 1 && ifocal >= 0, 'ifocal must be in range [0, 1]');
      assert(ifocal >= 0, 'ifocal must be >= 0');
      assert(isreal(kL), 'ratio of ifocal and ellipticity too small (try increasing ifocal)');
      assert(nscale > 0, 'Scale must be positive');
      
      if ellipticity < 1e-4 && strcmpi(parity, 'helical')
        warning(['Recommend using Bessel beam instead of ' ...
          'helical Mathieu beam with small ellipticity']);
      end
      
      % Calculate grid of points and use BscPmParaxial, this simulates
      % something a little closer to an SLM rather than directly PM
      % onto the hemisphere.
      
      % Aperture coordinates
      x = linspace(-1, 1, p.Results.Nres);
      y = linspace(-1, 1, p.Results.Nres);
      [X, Y] = meshgrid(x, y);
      R = sqrt(X.^2 + Y.^2);
      R(R == 0) = 1;
      dW = 1.0;
      
      % Scale by different amounts depending on mapping
      switch lower(p.Results.projection)
        case 'tanasin'
          Rd = tan(asin(R));
          X = X .* Rd ./ R;
          Y = Y .* Rd ./ R;
%           dW = cos(asin(R));
          
        case 'linear'
          %theta = atan(sqrt(x.^2 + y.^2));
          % Nothing to do
%           dW = 1.0;
          
        case 'tan'
          Rd = tan(R);
          X = X .* Rd ./ R;
          Y = Y .* Rd ./ R;
%           dW = 1./cos(atan(R)).^2;
          
        otherwise
          error('ott:BscPmParaxial:mapping', 'Unrecognized mapping value');
      end
      
      mask = R > 1;
      
      % Scale coordinates (related to theta_0)
      X = nscale .* X;
      Y = nscale .* Y;
      
      % Calculate coordinate range for eScalar calculation
      % TODO: Allow the user to set grid sizes
%       x = nscale*linspace(-1, 1, p.Results.Nres);
%       y = nscale*linspace(-1, 1, p.Results.Nres)';
      
      % TODO: Can we optimise the choice of coordinates or
      %     is Cartesian the best basis?
      %     Perhaps we could do it in elliptic and convert back to
      %     Cartesian instead, might be slightly more optimal?
      
      % Calculate elliptical coordinates
      [u,v]=BscPmMathieu.cart2ellip(ifocal,X,Y);
      
      % This is likely and can reduce problem size to 1/4
      % TODO: It might be better to change the coordinate system
      %     or implement the symmetry optimisations manually
      [uu, ~, icu] = unique(u);
      [uv, ~, icv] = unique(v);
      
      % Calculate scalar field
      switch parity
        case 'even'
          eU = BscPmMathieu.Mathieu_Je(uu, order, ellipticity);
          eV = BscPmMathieu.Mathieu_ce(uv, order, ellipticity);
          
          eScalar = reshape(eU(icu).*eV(icv), size(u));
        case 'odd'
          eU = BscPmMathieu.Mathieu_Jo(uu, order, ellipticity);
          eV = BscPmMathieu.Mathieu_se(uv, order, ellipticity);
                  
          eScalar = reshape(eU(icu).*eV(icv), size(u));
        case 'helical+'
          eU1 = BscPmMathieu.Mathieu_Je(uu, order, ellipticity);
          eV1 = BscPmMathieu.Mathieu_ce(uv, order, ellipticity);
          
          eU2 = BscPmMathieu.Mathieu_Jo(uu, order, ellipticity);
          eV2 = BscPmMathieu.Mathieu_se(uv, order, ellipticity);
          
          eScalar = reshape(eU1(icu).*eV1(icv) ...
              + 1i.*eU2(icu).*eV2(icv), size(u));
        case 'helical-'
          eU1 = BscPmMathieu.Mathieu_Je(uu, order, ellipticity);
          eV1 = BscPmMathieu.Mathieu_ce(uv, order, ellipticity);
          
          eU2 = BscPmMathieu.Mathieu_Jo(uu, order, ellipticity);
          eV2 = BscPmMathieu.Mathieu_se(uv, order, ellipticity);
          
          eScalar = reshape(eU1(icu).*eV1(icv) ...
              - 1i.*eU2(icu).*eV2(icv), size(u));
            
        otherwise
          error('Unsupported parity type');
      end
      
      % Scale power according to mapping function
      eScalar = eScalar .* dW;
      
      % Remove outside values
      eScalar(mask) = 0;
    end
  end
  
  methods (Static, Hidden)
    
    function G = Dmathieu_ce(x,r,q)
      %Angular Mathieu function type sin with an approximation of N=100
      %terms, it has three arguments, the first one is the variable, second
      %one is the ellipticity parameter q and the third one is the
      %eigenvalue m.
      %Note: Approximation can be improved by changing N

      N=100;
      G=0;
      Ar = BscPmMathieu.Fun_V_ce(r,q,N);
      
      
      if q==0
        
        G = -r*sin(r*x);
        
      else
        
        if r==0 || mod(r,2)==0
          
          for s=1:N
            
            G = G - Ar(s)*2*(s-1)*sin(2*(s-1)*x);
          end
          
        else
          
          for s=1:N
            
            G = G - Ar(s)*(2*(s-1)+1)*sin((2*(s-1)+1)*x);
          end
          
        end
        
      end
    end
    
    
    
    function G = Dmathieu_se(x,r,q)
      %Derivative of the angular sin type Mathieu function with an
      %approxiamtion of N=100 terms. It has 3 arguments, first one is the
      %variable, second one is the ellipticity parameter and the third one
      %is the eigenvalue.
      %Note: Approximation can be improved by changing N.
      
      N=100;
      G=0;
      B = BscPmMathieu.Fun_V_se(r,q,N);
      F1=0;
      x0 = 0.1;
      
      if q==0
        
        G = r*cos(r*x);
        
      else
        
        if mod(r,2)==0
          
          for s=1:N
            
            G = G + B(s)*2*s*cos(2*s*x);
            F1 = F1 + B(s)*sin(2*s*x0);
          end
          
        else
          
          for s=1:N
            
            G = G + B(s)*(2*(s-1)+1)*cos((2*(s-1)+1)*x);
            F1 = F1 + B(s)*sin((2*(s-1)+1)*x0);
          end
          
        end
        
        if F1<0
          G=-G;
        end
      end
    end
    
    
    
    
    function F = Dmathieu_Jo(x,r,q)
      %Radial Mathieu function type Se (or Jo) with an approximations of
      %N=100 terms. It has 3 arguments, first one is the
      %variable, second one is the ellipticity parameter and the third one
      %is the eigenvalue.
      %Note: Approximation can be improved by changing N.
      
      N=100;
      F=0;
      B = BscPmMathieu.Fun_V_se(r,q,N);
      arg = 2*sqrt(q)*cosh(x);
      
      
      if mod(r,2)==0
        
        Norm=BscPmMathieu.Dmathieu_se(pi/2,r,q)/(q*B(1));
        
        for s=1:N
          
          F = F + ((-1)^s)*2*s*B(s)*((((sech(x)).^2).*besselj(2*s,arg))+...
            sqrt(q)*tanh(x).*sinh(x).*(besselj(2*s-1,arg)-besselj(2*s+1,arg)));
        end
        
      else
        
        Norm=BscPmMathieu.Mathieu_se(pi/2,r,q)/(sqrt(q)*B(1));
        
        for s=1:N
          
          F = F + ((-1)^(s-1))*(2*(s-1)+1)*B(s)*((((sech(x)).^2).*besselj(2*(s-1)+1,arg))+...
            sqrt(q)*tanh(x).*sinh(x).*(besselj(2*(s-1),arg)-besselj(2*s,arg)));
        end
        
      end
      
      F = Norm*F;
    end

    function Gr = Dmathieu_Je(x,r,q)
      % Radial Mathieu function type Ce (or Je) with an approximation of
      % N=100 terms. It has 3 arguments, first one is the
      % variable, second one is the ellipticity parameter and the third one
      % is the eigenvalue.
      % Note: Approximation can be improved by changing N.
      
      N=100;
      F=0;
      Ar = BscPmMathieu.Fun_V_ce(r,q,N);
      arg = 2*sqrt(q)*cosh(x);
      
      
      if r==0 || mod(r,2)==0
        
        Norm=BscPmMathieu.Mathieu_ce(pi/2,r,q)/Ar(1);
        
        for s=1:N
          
          F = F + ((-1)^(s-1))*Ar(s)*(besselj(2*(s-1)-1,arg)-besselj(2*(s-1)+1,arg));
          
        end
        
      else
        
        Norm=BscPmMathieu.Dmathieu_ce(pi/2,r,q)/(sqrt(q)*Ar(1));
        
        for s=1:N
          
          F = F + ((-1)^s)*Ar(s)*(besselj(2*(s-1),arg)-besselj(2*s,arg));
        end
        
      end
      
      Gr = Norm*sqrt(q)*sinh(x).*F;
    end
    
    
    function M = Fun_MCnon(N,q)
      % Function that outputs the matrix M of odd coefficients for the calculation
      %of Mathieu functions type cosine. It carries 4 arguments; the first and the
      %second correspond to i and j respectively, the third argument is the
      %size of the matrix and the fourth argument is the parameter of ellipticity.

      for i=1:N    % rows
        for j=1:N % columns
          if j>=i
            if j==i
              if i==1
                M(i,j)=1+q;
              else
                M(i,j)=(2*(i-1) + 1)^2;
              end
            elseif j==i+1
              M(i,j)=q;
            else
              M(i,j)=0;
            end
          else
            M(i,j)=M(j,i);
          end
        end
      end
    end
    
    
    function M = Fun_MCpar(N,q)
      % Function that outputs the matrix M of even coefficients for the calculation
      %of Mathieu functions type cosine. It carries 4 arguments; the first and the
      %second correspond to i and j respectively, the third argument is the
      %size of the matrix and the fourth argument is the parameter of ellipticity.

      for i=1:N    % rows
        
        for j=1:N % columns
          
          if j>=i
            
            if j==i
              
              M(i,j)=(2*(i-1))^2;
              
            elseif j==i+1
              
              if i==1
                M(i,j)=sqrt(2)*q;
              else
                M(i,j)=q;
              end
              
            else
              M(i,j)=0;
            end
            
          else
            
            M(i,j)=M(j,i);
          end
        end
      end
    end

    function M = Fun_MSnon(N,q)
      % Function that outputs the matrix M of odd coefficients for the calculation
      %of Mathieu functions type sin. It carries 4 arguments; the first and the
      %second correspond to i and j respectively, the third argument is the
      %size of the matrix and the fourth argument is the parameter of ellipticity.

      for i=1:N    % rows
        for j=1:N % columns
          if j>=i
            if j==i
              if i==1
                M(i,j)=1-q;
              else
                M(i,j)=(2*(i-1) + 1)^2;
              end

            elseif j==i+1
              M(i,j)=q;

            else
              M(i,j)=0;
            end

          else
            M(i,j)=M(j,i);

          end
        end
      end
    end
    
    
    function M = Fun_MSpar(N,q)
      % Function that outputs the matrix M of even coefficients for the calculation
      %of Mathieu functions type sin. It carries 4 arguments; the first and the
      %second correspond to i and j respectively, the third argument is the
      %size of the matrix and the fourth argument is the parameter of ellipticity.

      for i=1:N    % rows
        
        for j=1:N % columns
          
          if j>=i
            
            if j==i
              
              M(i,j)=(2*i)^2;
              
            elseif j==i+1
              
              M(i,j)=q;
              
            else
              M(i,j)=0;
            end
            
          else
            
            M(i,j)=M(j,i);
          end
        end
      end
    end
    
    function  Ar = Fun_V_ce(r,q,N)
      % Calculate the k-th eigenvector of the angular Mathieu function
      % of type cosine with an approximation to order N

      if r==0
        M=BscPmMathieu.Fun_MCpar(N,q);
        r1=1;

      elseif   mod(r,2)==0
        M=BscPmMathieu.Fun_MCpar(N,q);
        r1=(r/2)+1;

      else
        M=BscPmMathieu.Fun_MCnon(N,q);
        r1=(r+1)/2;
      end

      [v,d]=eig(M);   % Calculate eigenvectors, v, and eigenvalues d
      [~,j]=sort(diag(d));  % Sort by eigenvalues, j contains new positions
      A = v(:, j);

      %D=d(r1);    % eigenvalue a_r
      A=A(:,r1);    % eigenvector

      if r==0 || mod(r,2)==0
        A(1)= A(1)/sqrt(2);
      end

      Ar=A;

      if Ar(1)<0
        Ar=-Ar;
      end
    end

    function  B = Fun_V_se(r,q,N)
      % Calculate the k-th eigenvector of the angular Mathieu function
      % of type sin with an approximation to order N.

      if mod(r,2)==0
        
        M=BscPmMathieu.Fun_MSpar(N,q);
        r1=r/2;
        
      else
        M=BscPmMathieu.Fun_MSnon(N,q);
        r1=(r+1)/2;
      end
      
      [v,d]=eig(M);   % Calculate eigenvectors, v, and eigenvalues, d
      [~,j]=sort(diag(d));  % Sort by eigenvalues, j contains new positions
      B = v(:, j);
      
      %D=d(r1);    % eigenvalue a_r
      B=B(:,r1);    % eigenvector
      
      if B(1)<0
        B=-B;
      end
    end

    function F = Mathieu_ce(x,r,q)
      %Angular Mathieu function type cosine with an approximations of N=100
      %terms, it has 3 arguments, the first one is the variable, the second
      %one is the eigenvalue and the third one is the ellipticity parameter.
      %Note: approximation can be improved by changing N.

      N=100;
      F=0;
      Ar = BscPmMathieu.Fun_V_ce(r,q,N);
      
      
      if r==0 || mod(r,2)==0
        
        for s=1:N
          
          F = F + Ar(s)*cos(2*(s-1)*x);
        end
        
      else
        
        for s=1:N
          
          F = F + Ar(s)*cos((2*(s-1)+1)*x);
        end
        
      end
    end
    
    
    function F = Mathieu_se(x,r,q)
      %Angular Mathieu function type sin with an approximations of N=100
      %terms, it has 3 arguments, the first one is the variable, the second
      %one is the eigenvalue and the third one is the ellipticity parameter.
      %Note: approximation can be improved by changing N.
      
      N=100;
      F=0;
      B = BscPmMathieu.Fun_V_se(r,q,N);
      F1=0;
      x0 = 0.1;
      
      if q==0
        
        F = sin(r*x);
        
      else
        
        if mod(r,2)==0
          
          for s=1:N
            
            F = F + B(s)*sin(2*s*x);
            F1 = F1 + B(s)*sin(2*s*x0);
            
          end
          
        else
          
          for s=1:N
            
            F = F + B(s)*sin((2*(s-1)+1)*x);
            F1 = F1 + B(s)*sin((2*(s-1)+1)*x0);
            
          end
          
        end
        
        if F1<0
          F=-F;
        end
      end
    end
    
    function F = Mathieu_Jo(x,r,q)
      %Radial Mathieu function type Se (or Jo) with an approximations of N=100
      %terms, it has 3 arguments, the first one is the variable, the second
      %one is the eigenvalue and the third one is the ellipticity parameter.
      %Note: approximation can be improved by changing N.

      N=100;
      F=0;
      B = BscPmMathieu.Fun_V_se(r,q,N);
      arg = 2*sqrt(q)*cosh(x);

      if mod(r,2)==0
        Norm=BscPmMathieu.Dmathieu_se(pi/2,r,q)/(q*B(1));
        for s=1:N
          F = F + ((-1)^s)*2*s*B(s)*besselj(2*s,arg);
        end
      else
        Norm=BscPmMathieu.Mathieu_se(pi/2,r,q)/(sqrt(q)*B(1));
        
        for s=1:N
          F = F + ((-1)^(s-1))*(2*(s-1)+1)*B(s)*besselj(2*(s-1)+1,arg);
        end
        
      end
      
      F = Norm*tanh(x).*F;
      
    end
    
    function Fr = Mathieu_Je(x,r,q)
      %Radial Mathieu function type Ce (or Je) with an approximations of N=100
      %terms, it has 3 arguments, the first one is the variable, the second
      %one is the eigenvalue and the third one is the ellipticity parameter.
      %Note: approximation can be improved by changing N.

      N=100;
      F=0;
      Ar = BscPmMathieu.Fun_V_ce(r,q,N);
      arg = 2*sqrt(q)*cosh(x);

      if r==0 || mod(r,2)==0
        
        Norm=BscPmMathieu.Mathieu_ce(pi/2,r,q)/Ar(1);
        
        for s=1:N
          F = F + ((-1)^(s-1))*Ar(s)*besselj(2*(s-1),arg);
        end
        
      else
        
        Norm=BscPmMathieu.Dmathieu_ce(pi/2,r,q)/(sqrt(q)*Ar(1));
        
        for s=1:N
          F = F + ((-1)^s)*Ar(s)*besselj(2*(s-1)+1,arg);
        end
        
      end
      
      Fr = Norm*F;
    end
  end

  methods
    function beam = BscPmMathieu(NA, order, parity, varargin)
      % Construct a new Mathieu beam via point-matching
      %
      % Usage:
      %   beam = BscPmMathieu(NA, order, parity, ...) constructs a new
      %   Mathieu beam with the default NA and scale.
      %
      % Parameters:
      %   - NA    -- (numeric) Numerical aperture of beam
      %   - order -- (numeric) Mathieu beam order
      %   - parity -- (enum) Mathieu beam parity (odd, even or helical+/-)
      %
      % Optional named arguments:
      %   - scale -- (numeric) Scaling factor for pattern.
      %     Default: 1.0
      %
      %   - ellipticity -- (numeric) Ellipticity of beam (q).
      %     This is related to the interfocal separation and the transverse
      %     spatial frequency by `q = f^2 k_T^2 / 4` where f is the semi-focal
      %     distance and k_T is the transverse spatial frequency.
      %     Default: 1.0
      %
      %   - ifocal -- (numeric) Interfocal separation.
      %     Default: 1.0
      %
      %   - Nres -- Resolution for paraxial field calculation.
      %     Default: 256.
      %
      % For additional named arguments see ott.BscPmParaxial.

      % Parse inputs
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('scale', 1.0);
      p.addParameter('ellipticity', 1.0);
      p.addParameter('ifocal', 1.0);
      p.addParameter('Nres', 256);
      p.addParameter('projection', 'linear');
      p.parse(varargin{:});

      % Check parity and order
      assert(isnumeric(order) && numel(order) == 1 && order == round(order), ...
          'order must be numeric integer');
      assert(any(strcmpi(parity, {'even', 'odd','helical+', 'helical-'})),...
          'parity must be even, odd or helical only');
      assert(p.Results.ellipticity >= 0, 'Ellipticity should be positive');

      % Calculate scalar field
      eScalar = BscPmMathieu.calculate_scalar_field(...
        parity, order, p.Results.scale, p.Results.ellipticity, ...
        p.Results.ifocal, 'Nres', p.Results.Nres, 'projection', p.Results.projection);

      % Initialise base class
      unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];
      beam = beam@ott.BscPmParaxial(NA, eScalar, unmatched{:});

      % Store properties
      beam.parity = parity;
      beam.order = order;
      beam.scale = p.Results.scale;
      beam.ellipticity = p.Results.ellipticity;
      beam.ifocal = p.Results.ifocal;
    end

  end
end

