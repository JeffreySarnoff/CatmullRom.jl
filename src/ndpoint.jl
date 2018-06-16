const coord_symbols = (Symbol.(map(x->string(x,"c"),split(lcletters, "")))...,);
# (:ac, :bc, :cc, :dc, :ec, :fc, :gc, :hc, ..., :pc, :qc, :rc, :sc, :tc, :uc, :vc, :wc, :xc, :yc, :zc)

const coordsyms = [
    (:xc,),
    (:xc, :yc),
    (:xc, :yc, :zc),
    (:xc, :yc, :zc, :wc),
    (:xc, :yc, :zc, :wc, :vc),
    (:xc, :yc, :zc, :wc, :vc, :uc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic, :hc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic, :hc, :gc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic, :hc, :gc, :fc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic, :hc, :gc, :fc, :ec),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic, :hc, :gc, :fc, :ec, :dc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic, :hc, :gc, :fc, :ec, :dc, :cc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic, :hc, :gc, :fc, :ec, :dc, :cc, :bc),
    (:xc, :yc, :zc, :wc, :vc, :uc, :tc, :sc, :rc, :qc, :pc,:oc, :nc, :mc, :lc, :kc, :jc, :ic, :hc, :gc, :fc, :ec, :dc, :cc, :bc, :ac),
    ];
    
    

    
