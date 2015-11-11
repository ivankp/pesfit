// Developed by Ivan Pogrebnyak, MSU
#include "pesfit.hh"

template<typename T> inline T sq(T x) noexcept { return x*x; }

namespace std {
  template <typename T1, typename T2>
  istream& operator>>(istream& in, pair<T1,T2>& p) {
    string s;
    in >> s;
    size_t sep = s.find(':');
    if (sep==string::npos) throw invalid_argument(
      cat('\"',s,"\": pair values must be delimited by \':\'"));
    stringstream (s.substr(0,sep)) >> p.first;
    stringstream (s.substr(sep+1)) >> p.second;
    return in;
  }
}

// options ------------------
senum(Fit,(none)(gaus)(cb))
Fit::type fit_;
senum(Out,(pdf)(root))
Out::type out_;

string ofname, cfname, wfname;
vector<string> ifname;
vector<Color_t> colors;
Int_t nbins;
pair<double,double> xrange;
bool logy, fix_alpha;
int prec;
vector<pair<string,pair<double,double>>> new_ws_ranges;
// --------------------------

// global -------------------
TFile *ofile;
workspace *ws;
seqmap<seqmap<val_err<double>>> stats;
TTree *tree;
TCanvas *canv;
TLatex *lbl;
// --------------------------

inline Double_t get_FWHM(const TH1* hist) noexcept {
  const Double_t half_max = hist->GetMaximum()/2;
  return hist->GetBinCenter(hist->FindLastBinAbove(half_max))
       - hist->GetBinCenter(hist->FindFirstBinAbove(half_max));
}
inline Double_t get_FWHM(const TGraph* gr) noexcept {
  const Double_t half_max = max(gr).second/2;
  return rfindx(gr,half_max) - lfindx(gr,half_max);
}

void make_hist(TH1*& hist, const char* name, const string& proc,
               double scale, const char* branch) {
  if (!tree->GetListOfBranches()->Contains(branch)) return;

  const string cmd1(cat(
    branch,"/1000>>hist(",nbins,",",xrange.first,",",xrange.second,")"));
  const string cmd2(
    "HGamEventInfoAuxDyn.crossSectionBRfilterEff"
    "*HGamEventInfoAuxDyn.weight"
    "*(HGamEventInfoAuxDyn.isPassed==1)");
  tree->Draw(cmd1.c_str(),cmd2.c_str());
  cout << endl << name
       << endl << cmd1
       << endl << cmd2 << endl;
  TH1 *temp = get<TH1>(gDirectory,"hist");
  temp->Scale(1000.*scale/temp->GetBinWidth(1));

  stats[name]["xsec_"+proc] = temp->Integral(0,temp->GetNbinsX()+1,"width");

  if (!hist)
    (hist = (TH1*)temp->Clone(name))->SetDirectory(0);
  else hist->Add(temp);
}

vector<FitResult> fit(const initializer_list<TH1*>& hs) {

  vector<FitResult> res;
  if (fit_==Fit::cb) res.reserve(hs.size());

  int i=0;
  Color_t color;
  for (auto it=hs.begin(), end=hs.end(); it!=end; ++it) {
    TH1 *hist = *it;
    if (!hist) continue;

    const char *name = hist->GetName();
    hist->SetStats(false);
    hist->SetLineWidth(2);
    color = colors[(i++) % colors.size()];
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);

    hist->SetXTitle("m_{#gamma#gamma} [GeV]");
    hist->SetYTitle("d#sigma/dm_{#gamma#gamma} [fb/GeV]");
    hist->SetTitleOffset(1.3,"Y");
    switch (out_) {
      case Out::pdf:
        if (it==hs.begin()) {
          hist->Draw();
          lbl->DrawLatex(.24,0.88-0.04*i,"Entries");
          lbl->DrawLatex(.32,0.88-0.04*i,"mean");
          lbl->DrawLatex(.40,0.88-0.04*i,"stdev");
        } else hist->Draw("same");
        break;
      case Out::root:
        hist->SetDirectory(ofile);
        break;
    }

    Double_t mean  = hist->GetMean(),
             stdev = hist->GetStdDev();

    auto &hstat = stats[name];
    hstat["hist_N"] = hist->GetEntries();
    hstat["hist_mean"] = {mean,hist->GetMeanError()};
    hstat["hist_stdev"] = {stdev,hist->GetStdDevError()};

    hstat["hist_window_mean"] = window_mean(hist,120,130);

    switch (fit_) { // FITTING +++++++++++++++++++++++++++++++++++++++
      case Fit::none: break;

      case Fit::gaus: {
        auto fit_res = hist->Fit("gaus","S");
        mean  = fit_res->Value(1);
        stdev = fit_res->Value(2);

        auto &hstat = stats[name];
        hstat["gaus_mean"] = {mean,fit_res->Error(1)};
        hstat["gaus_stdev"] = {stdev,fit_res->Error(2)};
      break; }

      case Fit::cb: {
        cout << "\033[32mFitting " << name << "\033[0m" << endl;
        auto fit_res = ws->fit(hist);
        switch (out_) {
          case Out::pdf:
            fit_res.second->Draw("same");
            break;
          case Out::root:
            auto *fit_gr = new TGraph(*fit_res.second);
            fit_gr->SetName(cat(hist->GetName(),"_fit").c_str());
            fit_gr->SetTitle(fit_gr->GetName());
            ofile->Add(fit_gr);
            break;
        }

        auto &hstat = stats[name];
        for (const string& varname : {
          "crys_alpha_bin0", "crys_norm_bin0", "fcb_bin0", "gaus_kappa_bin0",
          "gaus_mean_offset_bin0", "mean_offset_bin0", "sigma_offset_bin0"
        }) {
          auto *var = static_cast<RooRealVar*>(
            fit_res.first->floatParsFinal().find(varname.c_str()));
          hstat[varname] = {var->getVal(),var->getError()};
        }
        hstat["FWHM"] = get_FWHM(fit_res.second);

        if (fix_alpha) if (!strcmp(name,"nominal")) {
          auto *alpha = (*ws)->var("crys_alpha_bin0");
          alpha->setRange(alpha->getVal(),alpha->getVal());
        }

      res.emplace_back(move(fit_res.first));

      break; }
    }

    if (out_==Out::pdf) {
      auto lblp = lbl->DrawLatex(.12,0.84-0.04*i,hist->GetName());
      lblp->SetTextColor(color);
      lblp->DrawLatex(.24,0.84-0.04*i,cat(hist->GetEntries()).c_str());
      lblp->DrawLatex(.32,0.84-0.04*i,
        cat(fixed,setprecision(2),mean).c_str());
      lblp->DrawLatex(.40,0.84-0.04*i,
        cat(fixed,setprecision(2),stdev).c_str());
    }
  }
  if (out_==Out::pdf) if (i) canv->SaveAs(ofname.c_str());

  return res;
}

int main(int argc, char** argv)
{
  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("input,i", po::value(&ifname)->multitoken()->required(),
       "*input root file names")
      ("output,o", po::value(&ofname)->required(),
       "*output pdf or root file name")
      ("config,c", po::value(&cfname),
       "configuration file name")

      ("logy,l", po::bool_switch(&logy),
       "logarithmic Y axis")
      ("fit,f", po::value(&fit_)->default_value(Fit::none),
       cat("fit type: ",Fit::_str_all()).c_str())
      ("workspace,w", po::value(&wfname)->default_value("data/ws.root"),
       "ROOT file with RooWorkspace for CB fits")
      ("fix-alpha", po::bool_switch(&fix_alpha),
       "fix crys_alpha_bin0 parameter after nominal fit")
      ("xrange,x", po::value(&xrange)->default_value({105,140},"105:140"),
       "histograms\' X range")
      ("nbins,n", po::value(&nbins)->default_value(100),
       "histograms\' number of bins")
      ("prec", po::value(&prec)->default_value(-1),
       "summary table precision, -1 prints uncertainty")
      ("colors", po::value(&colors)->multitoken()->
        default_value(decltype(colors)({602,46}), "{602,46}"),
       "histograms\' colors")

      ("ws-setRange", po::value(&new_ws_ranges),
       "call RooWorkspace::setRange()")
    ;

    po::positional_options_description pos;
    pos.add("input",-1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(desc).positional(pos).run(), vm);
    if (argc == 1) {
      cout << desc << endl;
      return 0;
    }
    if (vm.count("config")) {
      po::store( po::parse_config_file<char>(
        vm["config"].as<string>().c_str(), desc), vm);
    }
    po::notify(vm);

    const string ofext = ofname.substr(ofname.rfind('.')+1);
    if (ofext=="pdf") out_ = Out::pdf;
    else if (ofext=="root") out_ = Out::root;
    else throw runtime_error(
      "Output file extension "+ofext+" is not pdf or root"
    );

  } catch (std::exception& e) {
    cerr << "\033[31mArgs: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // end options ---------------------------------------------------

  if (out_==Out::root) ofile = new TFile(ofname.c_str(),"recreate");

  ws = new workspace(wfname);
  for (const auto& range : new_ws_ranges)
    ws->setRange(range.first.c_str(),range.second.first,range.second.second);

  TH1 *nom=nullptr,
      *scale_down=nullptr, *scale_up=nullptr,
      *res_down=nullptr,   *res_up=nullptr;

  // LOOP over input files
  for (const string& f : ifname) {
    TFile *file = new TFile(f.c_str(),"read");
    if (file->IsZombie()) return 1;
    cout << "Data file: " << f << endl;
    tree = get<TTree>(file,"CollectionTree");

    const size_t slash = f.rfind('/')+1;
    const double xsecscale = 1./get<TH1>(file,
      ("CutFlow_"+f.substr(slash,f.find('.')-slash)+"_weighted").c_str()
    )->GetBinContent(3);

    // Regex for process identification
    static regex proc_re(".*[\\._]?(gg.|VBF|ttH|WH|ZH)[0-9]*[\\._]?.*",
                         regex_icase);
    smatch proc_match;
    if (!regex_match(f, proc_match, proc_re))
      throw runtime_error(cat("Filename \"",f,"\" does not specify process"));
    const string proc(proc_match.str(1));

    // protect from repeated processes
    static unordered_set<string> procs;
    if (!procs.emplace(proc).second)
      throw runtime_error(cat("File \"",f,"\" repeats process ",proc));

    // Make or add histograms
    make_hist(nom,"nominal",proc,xsecscale,
              "HGamEventInfoAuxDyn.m_yy");
    make_hist(scale_down,"scale_down",proc,xsecscale,
              "HGamEventInfo_EG_SCALE_ALL__1downAuxDyn.m_yy");
    make_hist(scale_up,"scale_up",proc,xsecscale,
              "HGamEventInfo_EG_SCALE_ALL__1upAuxDyn.m_yy");
    make_hist(res_down,"res_down",proc,xsecscale,
              "HGamEventInfo_EG_RESOLUTION_ALL__1downAuxDyn.m_yy");
    make_hist(res_up,"res_up",proc,xsecscale,
              "HGamEventInfo_EG_RESOLUTION_ALL__1upAuxDyn.m_yy");

    delete file;
  }
  cout << endl;

  // ---------------------------------------

  if (out_==Out::pdf) {
    canv = new TCanvas();
    canv->SetMargin(0.1,0.04,0.1,0.1);
    if (logy) canv->SetLogy();
    canv->SaveAs((ofname+'[').c_str());

    lbl = new TLatex();
    lbl->SetTextFont(43);
    lbl->SetTextSize(15);
    lbl->SetNDC();
  }

  TH2 *corr = fit({nom}).front()
    ->correlationHist("corr_mat");
  corr->SetTitle("Nominal signal fit correlation matrix");

  switch (out_) {
    case Out::pdf:
      if (logy) canv->SetLogy(false);
      gStyle->SetPaintTextFormat(".3f");
      canv->SetMargin(0.17,0.12,0.1,0.1);
      corr->SetStats(false);
      corr->SetMarkerSize(1.8);
      corr->Draw("COLZ TEXT");
      canv->SaveAs(ofname.c_str());
      canv->SetMargin(0.1,0.04,0.1,0.1);
      if (logy) canv->SetLogy(true);
      break;
    case Out::root:
      corr->SetDirectory(ofile);
      break;
  }

  fit({scale_down,scale_up});
  fit({res_down,res_up});

  // Print summary page *********************************************
  if (out_==Out::pdf) {
    canv->Clear();
    vector<unique_ptr<TPaveText>> txt;
    int i, n=6;
    txt.reserve(n);
    for (i=0; i<n; ++i) {
      TPaveText *pt;
      txt.emplace_back(pt = new TPaveText(float(i)/n,0.,float(i+1)/n,1.,"NBNDC"));
      pt->SetFillColor(0);
    }
    i=1;
    int m = 1;
    txt[0]->AddText("");
    for (auto& hist : stats) {
      txt[i]->AddText(hist.first.c_str());
      for (auto& var : hist.second) {
        if (i==1) {
          txt[0]->AddText(var.first.c_str());
          ++m;
        }
        stringstream ss;
        if (var.first.substr(0,var.first.find('_'))=="xsec") {
          ss << setprecision(3) << var.second.val;
        } else if (var.first.substr(var.first.rfind('_')+1)=="N") {
          ss << fixed << setprecision(0) << var.second.val;
        } else {
          if (prec==-1) var.second.print(ss," #pm ");
          else ss << scientific << setprecision(prec) << var.second.val;
        }
        txt[i]->AddText(ss.str().c_str());
      }
      ++i;
    }
    TLine *line = new TLine();
    for (i=0; i<n; ++i) {
      txt[i]->Draw();
      if (i) line->DrawLineNDC(float(i)/n,0.,float(i)/n,1.);
    }
    for (i=1; i<m; ++i) line->DrawLineNDC(0.,float(i)/m,1.,float(i)/m);
    canv->SaveAs(ofname.c_str());
  }
  // ****************************************************************

  if (fit_==Fit::cb && out_==Out::pdf) {
    double scale      = stats["nominal"   ]["mean_offset_bin0"].val;
    double scale_down = stats["scale_down"]["mean_offset_bin0"].val;
    double scale_up   = stats["scale_up"  ]["mean_offset_bin0"].val;
           scale_down = scale_down - scale;
           scale_up   = scale_up   - scale;
    double scale_sym  = (scale_up-scale_down)/2;

    double win        = stats["nominal"   ]["hist_window_mean"].val;
    double win_down   = stats["scale_down"]["hist_window_mean"].val;
    double win_up     = stats["scale_up"  ]["hist_window_mean"].val;
           win_down   = win_down - win;
           win_up     = win_up   - win;
    double win_sym    = (win_up-win_down)/2;

    double res        = stats["nominal"   ]["sigma_offset_bin0"].val;
    double res_down   = stats["res_down"  ]["sigma_offset_bin0"].val;
    double res_up     = stats["res_up"    ]["sigma_offset_bin0"].val;
           res_down   = res_down - res;
           res_up     = res_up   - res;
    double res_sym    = (res_up-res_down)/2;

    double fwhm       = stats["nominal"   ]["FWHM"].val;
    double fwhm_down  = stats["res_down"  ]["FWHM"].val;
    double fwhm_up    = stats["res_up"    ]["FWHM"].val;
           fwhm_down  = fwhm_down - fwhm;
           fwhm_up    = fwhm_up   - fwhm;
    double fwhm_sym   = (fwhm_up-fwhm_down)/2;

    canv->Clear();
    vector<unique_ptr<TPaveText>> txt;
    txt.reserve(6);
    for (int i=0; i<6; ++i) {
      TPaveText *pt;
      txt.emplace_back(pt = new TPaveText(i/6.,0.,(i+1)/6.,1.,"NBNDC"));
      pt->SetFillColor(0);
    }

    txt[0]->AddText("[GeV]");
    txt[1]->AddText("Scale");
    txt[2]->AddText("Window");
    txt[3]->AddText("Resolution");
    txt[4]->AddText("HWHM");
    txt[5]->AddText("FWHM/FWHM_{nom}");

    txt[0]->AddText("Nominal");
    txt[1]->AddText(Form("%.3f",scale));
    txt[2]->AddText(Form("%.3f",win));
    txt[3]->AddText(Form("%.3f",res));
    txt[4]->AddText(Form("%.3f",fwhm/2));
    txt[5]->AddText("1");

    txt[0]->AddText("Variation");
    txt[1]->AddText(Form("%.3f, +%.3f",scale_down,scale_up));
    txt[2]->AddText(Form("%.3f, +%.3f",win_down,win_up));
    txt[3]->AddText(Form("%.3f, +%.3f",res_down,res_up));
    txt[4]->AddText(Form("%.3f, +%.3f",fwhm_down/2,fwhm_up/2));
    txt[5]->AddText(Form("%.3f, +%.3f",fwhm_down/fwhm,fwhm_up/fwhm));

    txt[0]->AddText("Average Variation");
    txt[1]->AddText(Form("#pm %.3f",scale_sym));
    txt[2]->AddText(Form("#pm %.3f",win_sym));
    txt[3]->AddText(Form("#pm %.3f",res_sym));
    txt[4]->AddText(Form("#pm %.3f",fwhm_sym/2));
    txt[5]->AddText(Form("#pm %.3f",fwhm_sym/fwhm));

    TLine *line = new TLine();

    for (int i=0; i<6; ++i) {
      txt[i]->Draw();
      if (i) line->DrawLineNDC(i/6.,0.,i/6.,1.);
    }
    for (int i=1; i<4; ++i) line->DrawLineNDC(0.,i/4.,1.,i/4.);
    canv->SaveAs(ofname.c_str());
  }

  if (out_==Out::root) {
    ofile->cd();
    TTree *tree = new TTree("stats","stats");

    structmap(val_err<double>,hist_t,
      (nominal)(scale_down)(scale_up)(res_down)(res_up));

    seqmap<hist_t> tstats;
    for (auto& hist : stats)
      for (auto& var : hist.second)
        tstats[var.first][hist.first] = var.second;

    for (auto& stat : tstats)
      tree->Branch(stat.first.c_str(), &stat.second,
        "nominal[2]/D:scale_down[2]/D:scale_up[2]/D:res_down[2]/D:res_up[2]/D");

    tree->Fill();
  }

  switch (out_) {
    case Out::pdf:
      canv->SaveAs((ofname+']').c_str());
      delete canv;
      delete lbl;
      break;
    case Out::root:
      ofile->Write(0,TObject::kOverwrite);
      ofile->Close();
      delete ofile;
      break;
  }

  return 0;
}
