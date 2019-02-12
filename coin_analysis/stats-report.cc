#include <cstddef>
#include <vector>
#include <string>
#include <stdexcept>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TDatime.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TStyle.h>

#include <iostream>

#include "pg.hh"
#define START_DATE "2017-04-22"

struct points_t {
  std::vector<double> ts, vs;
};

struct points_t_error {
  std::vector<double> ts, vs, es;
};

points_t to_points(pg::untyped_result_set res)
{
  points_t pts;
  pts.ts.reserve(res.size()), pts.vs.reserve(res.size());
  for (size_t i=0; i<res.size(); ++i) {
    pts.ts.push_back(stod(res[i]["time"]));
    pts.vs.push_back(stod(res[i]["value"]));
  }
  return pts;
}

points_t_error to_points_error(pg::untyped_result_set res)
{
  points_t_error pts;
  pts.ts.reserve(res.size()), pts.vs.reserve(res.size()), pts.es.reserve(res.size());
  for (size_t i=0; i<res.size(); ++i) {
    pts.ts.push_back(stod(res[i]["time"]));
    pts.vs.push_back(stod(res[i]["value"]));
    pts.es.push_back(stod(res[i]["error"]));
  }
  return pts;
}

void draw(const std::string& title, const points_t& pts)
{
  TGraph g(pts.ts.size(), pts.ts.data(), pts.vs.data());
  g.SetTitle(title.c_str());
  gStyle->SetTitleFont(62);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetPadLeftMargin(0.06);
  gStyle->SetPadRightMargin(0.04);
  if(title=="el") {
    g.GetYaxis()->SetTitle("ns");
    g.GetYaxis()->SetRangeUser(0,3000000);
  }
  if (title=="efficiency"){
    g.GetYaxis()->SetRangeUser(0,1.22);
  }
  if (title=="AvgnS1"){
    g.GetYaxis()->SetRangeUser(0,10);
  }
  if (title=="AvgnS2"){
    g.GetYaxis()->SetRangeUser(0,10);
  }
  if (title=="noisetriggerrate"){
    g.GetYaxis()->SetRangeUser(0,1);
  }
  if (title=="sparkingrate"){
    g.GetYaxis()->SetTitle("Rate[Hz]");
    g.GetYaxis()->SetRangeUser(0,0.2);
  }
  g.SetMarkerColor(kRed);
  g.SetMarkerStyle(20);
  g.GetXaxis()->SetLabelOffset(0.03);
  g.GetXaxis()->SetTimeDisplay(1);
  g.GetXaxis()->SetTimeFormat("#splitline{%b-%d}{%Y}");
  g.GetXaxis()->SetTimeOffset(0,"cst");
  g.GetXaxis()->SetLabelSize(0.03);
  //TCanvas c("c","c",800,600);
  TCanvas c("c","c",1136,640);
  c.cd();
  g.Draw();
  c.Print((title+".png").c_str());
}

void draw_error(const std::string& title, const points_t_error& pts)
{
  TGraphErrors g(pts.ts.size(), pts.ts.data(), pts.vs.data(),0,pts.es.data());
  g.SetTitle(title.c_str());
  gStyle->SetTitleFont(62);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetPadLeftMargin(0.06);
  gStyle->SetPadRightMargin(0.04);
  g.SetMarkerColor(kRed);
  g.SetMarkerStyle(20);
  g.GetYaxis()->SetTitle("Rate[Hz]");
  g.GetXaxis()->SetLabelOffset(0.03);
  g.GetXaxis()->SetTimeDisplay(1);
  g.GetXaxis()->SetTimeFormat("#splitline{%b-%d}{%Y}");
  g.GetXaxis()->SetTimeOffset(0,"cst");
  g.GetXaxis()->SetLabelSize(0.03);
  TCanvas c("c","c",1136,640);
  c.cd();
  TLatex la;
  if (title=="Po218"){
    g.GetYaxis()->SetRangeUser(0,0.05);
  }
  if (title=="Rn222"){
    g.GetYaxis()->SetRangeUser(0,0.05);
  }
  if (title=="peak_236"){
    g.GetYaxis()->SetTitle("uBq/kg");
    g.GetYaxis()->SetRangeUser(0,500);
  }
  if (title=="peak_400"){
    g.GetYaxis()->SetTitle("uBq/kg");
    g.GetYaxis()->SetRangeUser(0,500);
  }
  if (title=="Xe131"){
    g.GetYaxis()->SetTitle("uBq/kg");
    g.GetYaxis()->SetRangeUser(0,500);
  }
  if (title=="Bi214"){
    g.GetYaxis()->SetRangeUser(0,0.02);
  }
  if (title=="low_energy_er_counts"){
    g.SetMarkerSize(0.7);
    c.SetGridx();
    c.SetGridy();
    c.SetLogy();
  }
  if (title=="Kr"){
    g.GetYaxis()->SetRangeUser(0,0.00002);
  }
  if (title=="Xe131"){
    g.GetYaxis()->SetRangeUser(0,20000);
  }
  if (title=="ss2count"){
    g.GetYaxis()->SetRangeUser(0,0.07);
  }
  if (title=="ss1candidate"){
    g.GetYaxis()->SetRangeUser(0,5);
  }
  if (title=="ss1count"){
//    g.GetYaxis()->SetRangeUser(0,0.003);
  }
  g.Draw();
  if (title=="low_energy_er_counts"){
    TDatime td(pts.ts.data()[pts.ts.size()-1]);
    la.DrawLatex(1.472e9,1e-2,Form("Last Point: %d",td.GetDate()));
    la.DrawLatex(1.472e9,0.5e-2,Form("%f",pts.vs.data()[pts.ts.size()-1]));
  }
  c.Print((title+".png").c_str());
}

points_t reduce_by_avg(pg::connection_t con, const std::string& param, const std::string& window)
{
  auto res = pg::query(
    con,
    "select "
    "  extract(epoch from date_trunc_interval(file_open_time, $2::interval)) as time, "
    "  avg((params->$1)::double precision) as value "
    "from files, runs "
    "where files.run_no=runs.run_no "
    "      and ((params->$1)::double precision <> 'NaN') "
    "      and ((params->$1)::double precision <> 'Infinity') "
    "      and ((params->$1)::double precision <> '-Infinity') "
    "      and quality "
    "      and files.file_open_time > $3::timestamp "
    "      and type='PHYSICS_DM'"   
    "group by 1 "
    "order by time",
    param, window, START_DATE);

  return to_points(res);
}

points_t_error reduce_as_rate(pg::connection_t con, const std::string& param, const std::string& window)
{
  auto res = pg::query(
    con,
    "select "
    "  extract(epoch from date_trunc_interval(file_open_time, $2::interval)) as time, "
    "  sum((params->$1)::double precision)"
    "  / sum(extract(epoch from (file_close_time-file_open_time))) as value, "
    "  sqrt(greatest(sum((params->$1)::double precision),0))"
    "  / sum(extract(epoch from (file_close_time-file_open_time))) as error "
    "from files,runs "
    "where files.run_no=runs.run_no "
    "      and ((params->$1)::double precision <> 'NaN') "
    "      and ((params->$1)::double precision <> 'Infinity') "
    "      and ((params->$1)::double precision <> '-Infinity') "
    "      and quality "
    "      and files.file_open_time > $3::timestamp "
    "      and type='PHYSICS_DM'"   
    "group by 1 "
    "order by time",
    param, window, START_DATE);

  return to_points_error(res);
}

points_t plain_pts(pg::connection_t con, const std::string& param)
{
  auto res = pg::query(
    con,
    "select "
    "  extract(epoch from file_open_time) as time, "
    "  (params->$1)::double precision as value "
    "from files "
    "where params ? $1 "
    "      and ((params->$1)::double precision <> 'NaN') "
    "      and ((params->$1)::double precision <> 'Infinity') "
    "      and ((params->$1)::double precision <> '-Infinity') "
    "      and ((params->$1)::double precision <10000000)"
    "      and files.file_open_time > $2::timestamp "
    "      and quality "
    "order by time",
    param, START_DATE);

  return to_points(res);
}

int main(int argc, char *argv[])
{
  auto con = pg::connect(getenv("DAQ_DBURL"));

  if (!(argc==3 || argc==4)) {
    std::cerr << "Usage: " << argv[0] << " type param [window]\n";
    return 1;
  }

  std::string type = argv[1], param = argv[2], window;
  if (argc == 4) window = argv[3];

#define ARGASSERT(good, msg)                    \
  if (!(good)) {                                \
    std::cerr << msg << '\n';                   \
    return 1;                                   \
  }

  if (type == "rate") {
    ARGASSERT(argc==4, "Type rate requires window argument");
    draw_error(param, reduce_as_rate(con, param, window));
  } else if (type == "avg") {
    ARGASSERT(argc==4, "Type avg requires window argument");
    draw(param, reduce_by_avg(con, param, window));
  } else if (type == "plain") {
    ARGASSERT(argc==3, "Type plain requires no window argument");
    draw(param, plain_pts(con, param));
  } else {
    ARGASSERT(false, "Unknown type");
  }
  
  
  return 0;
}
