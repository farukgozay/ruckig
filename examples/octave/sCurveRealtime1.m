% sCurveRealtime1.m - Asymmetric 3-Phase S-Curve generator
% Produces acceleration → cruise → deceleration with separate a_max and d_max

function [pos, vel, acc, dec, jerk] = sCurveRealtime1(TargetPos, CurrentPos, Ts, v_max, a_max, d_max, jerk_max)
    % keep state between calls
    persistent D_rem T_brake D_brake prev_vel t_current q0 qf signD
    persistent D_acc D_cruise D_dec T_acc T_cruise T_dec tf
    persistent coeffsA coeffsD coeffsA1 coeffsA2 v_init a_init reset_flag
    persistent v_peak v_peak_ideal v0 a0

    if isempty(prev_vel)
        prev_vel = 0;
    end

    % first call or after reset, set up initial conditions
    if isempty(reset_flag)
        reset_flag = true;
        v_init     = 0;
        a_init     = 0;
        t_current  = 0;
    end

    if isempty(t_current)
        t_current = 0;
    end

    % when a new target arrives or profile is done, recompute everything
    if reset_flag || TargetPos ~= qf
        reset_flag = false;
        t_current  = 0;
        q0         = CurrentPos;
        qf         = TargetPos;
        signD      = sign(qf - q0);

        v0 = v_init;
        a0 = a_init;

        D = abs(qf - q0);
        v_peak_ideal = sqrt((2*D*a_max*d_max)/(a_max + d_max));

        if sign(v0) == signD && v0 == 0
            % standard single-phase accel
            if v_peak_ideal <= v_max
                v_peak   = v_peak_ideal;
                D_acc    = (v_peak^2 - v0^2)/(2*a_max);
                D_dec    = v_peak^2/(2*d_max);
                D_cruise = 0;
            else
                v_peak   = v_max;
                D_acc    = (v_peak^2 - v0^2)/(2*a_max);
                D_dec    = v_peak^2/(2*d_max);
                D_cruise = D - (D_acc + D_dec);
            end

            T_acc    = (v_peak - v0)/a_max;
            T_dec    = v_peak/d_max;
            T_cruise = max(0, D_cruise/v_peak);
            coeffsA  = quintic_coeffs(q0, v0, a0, ...
                                      q0 + signD*D_acc, signD*v_peak, 0, T_acc);
        else
            % interrupt: brake to zero, then accel from zero
            D_brake     = v0^2/(2*d_max);
            T_brake     = abs(v0)/d_max;
            D_rem       = D - D_brake;

            v_peak      = min(v_peak_ideal, v_max);
            D_accel     = v_peak^2/(2*a_max);
            T_accel     = v_peak/a_max;
            D_decel     = v_peak^2/(2*d_max);
            T_decel     = v_peak/d_max;

            D_cruise    = max(0, D_rem - (D_accel + D_decel));
            T_cruise    = D_cruise / v_peak;

            % ensure at least one time step of cruise
            T_cruise    = max(T_cruise, Ts);
            D_cruise    = v_peak * T_cruise;

            T_acc       = T_brake + T_accel;
            tf          = T_acc + T_cruise + T_decel;

            q1          = q0 + sign(v0)*D_brake;
            coeffsA1    = quintic_coeffs(q0, v0, a0, q1, 0, 0, T_brake);

            q2          = q1 + signD*D_accel;
            coeffsA2    = quintic_coeffs(q1, 0, 0, q2, signD*v_peak, 0, T_accel);
        end

        % end decel polynomial
        coeffsD = quintic_coeffs(q0 + signD*(D_acc + D_cruise), ...
                                 signD*v_peak, 0, qf, 0, 0, T_dec);
    end

    % compute current phase at t_current
    if sign(v0) == signD && v0 == 0
        if t_current <= T_acc
            t   = t_current;
            pos = polyval(flip(coeffsA'), t);
            vel = polyval(polyder(flip(coeffsA')), t);
            acc = polyval(polyder(polyder(flip(coeffsA'))), t);
            jerk = polyval(polyder(polyder(polyder(flip(coeffsA')))), t);
            dec = 0;
        elseif t_current <= T_acc + T_cruise
            t_rel = t_current - T_acc;
            pos   = q0 + signD*(D_acc + v_peak*t_rel);
            vel   = signD*v_peak;
            acc   = 0;
            jerk  = 0;
            dec   = 0;
        else
            t_rel = t_current - (T_acc + T_cruise);
            pos   = polyval(flip(coeffsD'), t_rel);
            vel   = polyval(polyder(flip(coeffsD')), t_rel);
            acc   = polyval(polyder(polyder(flip(coeffsD'))), t_rel);
            jerk  = polyval(polyder(polyder(polyder(flip(coeffsD')))), t_rel);
            dec   = max(0, -acc);
        end
    else
        if t_current <= T_brake
            t   = t_current;
            pos = polyval(flip(coeffsA1'), t);
            vel = polyval(polyder(flip(coeffsA1')), t);
            acc = polyval(polyder(polyder(flip(coeffsA1'))), t);
            jerk = polyval(polyder(polyder(polyder(flip(coeffsA1')))), t);
            dec = max(0, -acc);
        elseif t_current <= T_acc
            t   = t_current - T_brake;
            pos = polyval(flip(coeffsA2'), t);
            vel = polyval(polyder(flip(coeffsA2')), t);
            acc = polyval(polyder(polyder(flip(coeffsA2'))), t);
            jerk = polyval(polyder(polyder(polyder(flip(coeffsA2')))), t);
            dec = 0;
        elseif t_current <= T_acc + T_cruise
            t_rel = t_current - T_acc;
            pos   = q0 + signD*(D_acc + v_peak*t_rel);
            vel   = signD*v_peak;
            acc   = 0;
            jerk  = 0;
            dec   = 0;
        else
            t_rel = t_current - (T_acc + T_cruise);
            pos   = polyval(flip(coeffsD'), t_rel);
            vel   = polyval(polyder(flip(coeffsD')), t_rel);
            acc   = polyval(polyder(polyder(flip(coeffsD'))), t_rel);
            jerk  = polyval(polyder(polyder(polyder(flip(coeffsD')))), t_rel);
            dec   = max(0, -acc);
        end
    end

    % enforce limits
    vel  = max(min(vel,  v_peak),   -v_peak);
    acc  = max(min(acc,  a_max),    -a_max);
    jerk = max(min(jerk, jerk_max), -jerk_max);

    % advance time, reset when done
    t_current = t_current + Ts;
    if t_current > tf
        pos        = qf;
        vel        = 0;
        acc        = 0;
        dec        = 0;
        jerk       = 0;
        reset_flag = true;
        v_init     = vel;
        a_init     = acc;
    end

    prev_vel = vel;
end

function c = quintic_coeffs(q0, v0, a0, qf, vf, af, T)
    % solve 5th-degree polynomial coefficients
    A = [1,0,0,0,0,0;
         0,1,0,0,0,0;
         0,0,2,0,0,0;
         1,T,T^2,T^3,T^4,T^5;
         0,1,2*T,3*T^2,4*T^3,5*T^4;
         0,0,2,6*T,12*T^2,20*T^3];
    b = [q0; v0; a0; qf; vf; af];
    c = A\b;
end
