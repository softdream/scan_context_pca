#include "laserSimulation.h"

#include "scan_context_pca.h"

#include "scan_context.h"

#include "KDTreeVectorOfVectorsAdaptor.h"
#include <vector>

#include <opencv2/opencv.hpp>

static const unsigned char r[64] = { 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   10,  20,  30,  40,  50,  60,  70,  80,  90,  100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };
static const unsigned char g[64] = { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 240, 230, 220, 210, 200, 190, 180, 170, 160, 150, 140, 130, 120, 110, 100, 90,  80,  70,  60,  50,  40,  30,  20,  10,  0 };
static const unsigned char b[64] = { 255, 240, 220, 200, 180, 160, 140, 120, 100, 80,  60,  40,  20,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 };


template<int Dimension>
using myVectors = std::vector<Eigen::Matrix<float, Dimension, 1>>;

myVectors<20> vec;



        template<typename TT>
        const TT rad2deg( const TT radians)
        {
                return radians * 180.0 / M_PI;
        }

        template<typename TT>
        const TT deg2rad( const TT angle )
        {
                return angle * M_PI / 180.0;
        }

void drawABin(cv::Mat &image,  int ring_idx, const int sctor_idx)
{
        float angleStart1 = (sctor_idx) * 6 * M_PI / 180;
        float radiusStart1 = (ring_idx) * 20;
        float angleStart2 = (sctor_idx + 1) * 6 * M_PI / 180;
        float radiusStart2 = (ring_idx + 1) * 20;

        std::vector<std::vector<cv::Point>> ppPoints;
        std::vector<cv::Point> pPoints;

        for (int i = 0; i <= 6; i++) {
                cv::Point p( 450 - radiusStart1 * sin(angleStart1 + i * 0.017453293),
                             450 - radiusStart1 * cos(angleStart1 + i * 0.017453293));
                pPoints.push_back( p );
        }

        for (int i = 0; i <= 6; i++) {
                cv::Point p( 450 - radiusStart2 * sin( angleStart2 - i * 0.017453293),
                             450 - radiusStart2 * cos( angleStart2 - i * 0.017453293));
                pPoints.push_back( p );
        }

        ppPoints.push_back( pPoints );

        /*for (auto it : pPoints) {
                cv::circle(image, it, 3, cv::Scalar(0, 0, 255), -1);
                cv::imshow("scan distribution", image);
                cv::waitKey(2000);
        }*/

        cv::fillPoly( image, ppPoints, cv::Scalar( 0, 0, 255 ) );
        //cv::fillPoly( image, ppPoints, cv::Scalar(b[64 - count], g[64 - count], r[64 - count]) );
}


void displayScanDistribution( const sensor::LaserScan &scan, const int num = 0 )
{
        cv::Mat image = cv::Mat::zeros( 900, 900, CV_8UC3 );

        for( int i = 0; i < 20; i ++ ){
                cv::circle( image, cv::Point2f( 450, 450 ), 20 * ( i + 1 ), cv::Scalar( 0, 255, 0 ), 1 );
        }

        for( int i = 0; i < 60; i ++ ){
                float angle =  ( 6 * i );

                cv::Point2f endPoint( 450 - ::sin( angle * M_PI / 180 ) * 400,
                                      450 - ::cos( angle * M_PI / 180 ) * 400 );
                cv::line( image, cv::Point2f( 450, 450 ), endPoint, cv::Scalar( 0, 255, 0 ), 1 );
        }
        cv::arrowedLine( image, cv::Point2f( 450, 900 ), cv::Point2f( 450, 20 ), cv::Scalar( 255, 0, 0 ), 1 );
        cv::arrowedLine( image, cv::Point2f( 0, 450 ), cv::Point2f( 880, 450 ), cv::Scalar( 255, 0, 0 ), 1 );


        int ring_idx = 0, sctor_idx = 0;
        float radians = -3.14159f;

        for( int i = 0; i < scan.size(); i ++ ){
                float dist = scan.ranges[i];
                float angle = rad2deg<float>( radians ) + 180.0f;

                if( dist >= 0.009999998f && dist <= 10.0000000000f ){
                        ring_idx = std::max( std::min( 20 - 1, int(ceil( (dist / 10) * 20 )) ), 0 );

                        sctor_idx = std::max( std::min( 60 - 1, int(ceil( (angle / 360.0) * 60 )) ), 0 );

                        // draw the bin
                        //cv::Point2f binPoint( 450 - ( ring_idx * 20 ) * ::sin( sctor_idx * 6 * M_PI / 180 ) + 10, 
                        //                    450 - ( ring_idx * 20 ) * ::cos( sctor_idx * 6 * M_PI / 180 ) - 10 );
                        //cv::circle( image, binPoint, 3, cv::Scalar( 0, 0, 255 ), -1 );

                        drawABin( image, ring_idx, sctor_idx );
                }

                radians += 0.0043542264f;
        }


        cv::imshow( "scan_distribution_+" + std::to_string( num ), image );
}

void displayAScancontext( const Eigen::Matrix<float, 20, 60> &desc, const int num = 0)
{
        cv::Mat image = cv::Mat::zeros( 20 * 10, 60 * 10, CV_8UC3 );

        for( int i = 0; i < 20; i ++ ){
                cv::line( image, cv::Point( 0, 10 * i ), cv::Point( 60 * 10, 10 * i ), cv::Scalar( 67, 128, 94 ), 1 );
        }

        for( int i = 0; i < 60; i ++ ){
                cv::line( image, cv::Point( 10 * i, 0 ), cv::Point( 10 * i, 20 * 10 ), cv::Scalar( 67, 128, 94 ), 1 );
        }

        for( int i = 0; i < desc.rows(); i ++ ){ 
                for( int j = 0; j < desc.cols(); j ++ ){
                        if( desc( i, j ) > 0 ){
                                //std::cout<<"point: ( "<<i <<", "<<j <<" )"<<std::endl;
                                cv::Point points[1][4];
                                points[0][0] = cv::Point(j * 10, i * 10); 
                                points[0][1] = cv::Point(j * 10 + 10, i * 10); 
                                points[0][2] = cv::Point(j * 10 + 10, i * 10 + 10);
                                points[0][3] = cv::Point(j * 10, i * 10 + 10);   

                                const cv::Point* ppt[1] = { points[0] };
                                int npt[] = { 4 };

                                int count = static_cast<int>( desc( i, j ) );
                                cv::fillPoly(image, ppt, npt, 1, cv::Scalar(b[64 - count], g[64 - count], r[64 - count]));
                        }
                }
        }

        cv::imshow( "scancontext" + std::to_string( num ), image );

}


int main()
{
	simulation::Simulation simulation;
        simulation.openSimulationFile( "laser_data.txt" );
	
	int count = 0;
        while( count < 1450 ){
		sensor::LaserScan scan;
                simulation.readAFrameData( scan );
                std::cout<<"-------------- "<<simulation.getFrameCount()<<" -------------"<<std::endl;
	
		scancontext::ScanContextPCA<float> pca;
		Eigen::Matrix<float, 20, 1> key = pca.caculateResult( scan );
		std::cout<<"PCA Vector : "<<std::endl<<key<<std::endl;
		//context::ScanContext<float> sc;
		//Eigen::Matrix<float, 20, 1> key = sc.caculateResult( scan );
		//std::cout<<"Key Vector : "<<std::endl<<key<<std::endl;
		
		if( count == 550 ){
			displayAScancontext( pca.getDesc(), count );
			displayScanDistribution( scan, count );
		}
		
		if( count == 1071 ){
                        displayAScancontext( pca.getDesc(), count );
			displayScanDistribution( scan, count );
                }

		if( count == 1350 ){
                        displayAScancontext( pca.getDesc(), count );
			displayScanDistribution( scan, count );
			cv::waitKey(0);
                }

		//cv::waitKey(60);

		vec.push_back( key );	

		count ++;
	}

	// ----------------------------------------------- //
	Eigen::Matrix<float, 20, 1> query_pt;
	query_pt << 0,
         0,
  -2.78333,
     42.95,
  -2.36667,
     -1.15,
 -0.233333,
-0.0833333,
        -1,
 -0.216667,
 -0.116667,
      -0.1,
     -0.05,
         0,
         0,
         0,
         0,
         0,
         0,
         0;



	typedef KDTreeVectorOfVectorsAdaptor< myVectors<20>, float >  my_kd_tree_t;
        my_kd_tree_t mat_index( 20, vec, 10 );
        mat_index.index->buildIndex();

        // do a knn search
        const size_t num_results = 10;
        std::vector<size_t>   ret_indexes(num_results);
        std::vector<double> out_dists_sqr(num_results);

        nanoflann::KNNResultSet<double> resultSet(num_results);

        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

        std::cout << "knnSearch(nn=" << num_results << "): \n";
        for (size_t i = 0; i < num_results; i++)
                std::cout << "ret_index[" << i << "]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << std::endl;

	return 0;
}








