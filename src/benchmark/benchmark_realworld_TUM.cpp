#include "tools.hpp"
#include <ros/ros.h>
#include <Eigen/Eigenvalues>
#include <sensor_msgs/PointCloud2.h>
#include <pcl_conversions/pcl_conversions.h>
#include <geometry_msgs/PoseArray.h>
#include <random>
#include <ctime>
#include <tf/transform_broadcaster.h>
#include "bavoxel.hpp"

#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <malloc.h>

using namespace std;

template <typename T>
void pub_pl_func(T &pl, ros::Publisher &pub)
{
  pl.height = 1; pl.width = pl.size();
  sensor_msgs::PointCloud2 output;
  pcl::toROSMsg(pl, output);
  output.header.frame_id = "camera_init";
  output.header.stamp = ros::Time::now();
  pub.publish(output);
}

ros::Publisher pub_path, pub_test, pub_show, pub_cute;

int load_evo_pose_with_time(
    const std::string &pose_file,
    std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>> &pose_list,
    std::vector<double> &time_list) {
  time_list.clear();
  pose_list.clear();
  std::ifstream fin(pose_file);
  std::string line;
  Eigen::Matrix<double, 1, 7> temp_matrix;
  while (getline(fin, line)) {
    std::istringstream sin(line);
    std::vector<std::string> Waypoints;
    std::string info;
    int number = 0;
    while (getline(sin, info, ' ')) {
      if (number == 0) {
        double time;
        std::stringstream data;
        data << info;
        data >> time;
        time_list.push_back(time);
        number++;
      } else {
        double p;
        std::stringstream data;
        data << info;
        data >> p;
        temp_matrix[number - 1] = p;
        if (number == 7) {
          Eigen::Vector3d translation(temp_matrix[0], temp_matrix[1],
                                      temp_matrix[2]);
          Eigen::Quaterniond q(temp_matrix[6], temp_matrix[3], temp_matrix[4],
                               temp_matrix[5]);
          std::pair<Eigen::Vector3d, Eigen::Matrix3d> single_pose;
          single_pose.first = translation;
          single_pose.second = q.toRotationMatrix();
          pose_list.push_back(single_pose);
        }
        number++;
      }
    }
  }
  return pose_list.size();
}

void read_file(vector<IMUST> &x_buf, vector<pcl::PointCloud<PointType>::Ptr> &pl_fulls, string &prename)
{
  std::string pcds_dir = prename;
  std::string pose_file = prename + "/utm_L_opt_pose_W.txt";

  std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>> pose_list;
  std::vector<double> time_list;
  int pose_size = load_evo_pose_with_time(pose_file, pose_list, time_list);
  
  pcl::PCDReader reader;
  for(int m=0; m<pose_size; m++)
  {
    // Load point cloud from pcd file
    std::stringstream ss;
    ss << pcds_dir << "/" << std::setfill('0') << std::setw(6) << m << ".pcd";
    std::string filename = ss.str();

    pcl::PointCloud<PointType>::Ptr pl_ptr(new pcl::PointCloud<PointType>());
    pcl::PointCloud<pcl::PointXYZI> pl_tem;
    pcl::io::loadPCDFile(filename, pl_tem);
    for(pcl::PointXYZI &pp: pl_tem.points)
    {
      PointType ap;
      ap.x = pp.x; ap.y = pp.y; ap.z = pp.z;
      ap.intensity = pp.intensity;
      pl_ptr->push_back(ap);
    }

    pl_fulls.push_back(pl_ptr);

    IMUST curr;
    curr.R = pose_list[m].second; curr.p = pose_list[m].first; curr.t = time_list[m];
    x_buf.push_back(curr);

    if (x_buf.size() > 10)
      break;
  }
}

void data_show(vector<IMUST> x_buf, vector<pcl::PointCloud<PointType>::Ptr> &pl_fulls)
{
  IMUST es0 = x_buf[0];
  for(uint i=0; i<x_buf.size(); i++)
  {
    x_buf[i].p = es0.R.transpose() * (x_buf[i].p - es0.p);
    x_buf[i].R = es0.R.transpose() * x_buf[i].R;
  }

  pcl::PointCloud<PointType> pl_send, pl_path;
  int winsize = x_buf.size();
  for(int i=0; i<winsize; i++)
  {
    pcl::PointCloud<PointType> pl_tem = *pl_fulls[i];
    down_sampling_voxel(pl_tem, 0.05);
    pl_transform(pl_tem, x_buf[i]);
    pl_send += pl_tem;

    if((i%200==0 && i!=0) || i == winsize-1)
    {
      pub_pl_func(pl_send, pub_show);
      pl_send.clear();
      sleep(0.5);
    }

    PointType ap;
    ap.x = x_buf[i].p.x();
    ap.y = x_buf[i].p.y();
    ap.z = x_buf[i].p.z();
    ap.curvature = i;
    pl_path.push_back(ap);
  }

  pub_pl_func(pl_path, pub_path);
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "benchmark2");
  ros::NodeHandle n;
  pub_test = n.advertise<sensor_msgs::PointCloud2>("/map_test", 100);
  pub_path = n.advertise<sensor_msgs::PointCloud2>("/map_path", 100);
  pub_show = n.advertise<sensor_msgs::PointCloud2>("/map_show", 100);
  pub_cute = n.advertise<sensor_msgs::PointCloud2>("/map_cute", 100);

  // 读取数据
  string prename, ofname;
  vector<IMUST> x_buf;
  vector<pcl::PointCloud<PointType>::Ptr> pl_fulls;

  n.param<double>("voxel_size", voxel_size, 1);
  string file_path;
  n.param<string>("file_path", file_path, "");
  double pose_noise;
  // n.param<string>("pose_noise", pose_noise, "");

  read_file(x_buf, pl_fulls, file_path);

  // 处理数据
  IMUST es0 = x_buf[0];
  for(uint i=0; i<x_buf.size(); i++)
  {
    x_buf[i].p = es0.R.transpose() * (x_buf[i].p - es0.p);
    x_buf[i].R = es0.R.transpose() * x_buf[i].R;
  }

  win_size = x_buf.size();
  printf("The size of poses: %d\n", win_size);

  data_show(x_buf, pl_fulls);
  printf("Check the point cloud with the initial poses.\n");
  printf("If no problem, input '1' to continue or '0' to exit...\n");
  int a; cin >> a; if(a==0) exit(0);

  pcl::PointCloud<PointType> pl_full, pl_surf, pl_path, pl_send;
  for(int iterCount=0; iterCount<1; iterCount++)
  { 
    unordered_map<VOXEL_LOC, OCTO_TREE_ROOT*> surf_map;

    eigen_value_array[0] = 1.0 / 16;
    eigen_value_array[1] = 1.0 / 16;
    eigen_value_array[2] = 1.0 / 9;

    for(int i=0; i<win_size; i++)
      cut_voxel(surf_map, *pl_fulls[i], x_buf[i], i);

    pcl::PointCloud<PointType> pl_send;
    pub_pl_func(pl_send, pub_show);

    pcl::PointCloud<PointType> pl_cent; pl_send.clear();
    VOX_HESS voxhess;
    for(auto iter=surf_map.begin(); iter!=surf_map.end() && n.ok(); iter++)
    {
      iter->second->recut(win_size);
      iter->second->tras_opt(voxhess, win_size);
      iter->second->tras_display(pl_send, win_size);
    }

    pub_pl_func(pl_send, pub_cute);
    printf("\nThe planes (point association) cut by adaptive voxelization.\n");
    printf("If the planes are too few, the optimization will be degenerated and fail.\n");
    printf("If no problem, input '1' to continue or '0' to exit...\n");
    int a; cin >> a; if(a==0) exit(0);
    pl_send.clear(); pub_pl_func(pl_send, pub_cute);

    if(voxhess.plvec_voxels.size() < 3 * x_buf.size())
    {
      printf("Initial error too large.\n");
      printf("Please loose plane determination criteria for more planes.\n");
      printf("The optimization is terminated.\n");
      exit(0);
    }

    BALM2 opt_lsv;
    opt_lsv.damping_iter(x_buf, voxhess);

    for(auto iter=surf_map.begin(); iter!=surf_map.end();)
    {
      delete iter->second;
      surf_map.erase(iter++);
    }
    surf_map.clear();

    malloc_trim(0);
  }

  printf("\nRefined point cloud is publishing...\n");
  malloc_trim(0);
  data_show(x_buf, pl_fulls);
  printf("\nRefined point cloud is published.\n");

  ros::spin();
  return 0;

}


