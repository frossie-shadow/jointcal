// -*- C++ -*-
// 
#ifndef ASSOCIATIONS__H
#define ASSOCIATIONS__H

#include <string>
#include <iostream>

//#include "stringlist.h"

#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/geom/Box.h"

#include "lsst/meas/simastrom/RefStar.h"
#include "lsst/meas/simastrom/FittedStar.h"
#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/Point.h"

namespace lsst {
namespace meas {
namespace simastrom {

//! The class that implements the relations between MeasuredStar and FittedStar.
class Associations {
  public: 

  CcdImageList ccdImageList; // the catalog handlers
  RefStarList refStarList;// the (e.g.) USNO stars
  RefStarList photRefStarList; // the (e.g) Landolt stars
  FittedStarList fittedStarList; //  the std::list of stars that are going to be fitted

    // fit cuts and stuff:
  Point commonTangentPoint;

//  const CatalogLoader * load_it;

  void AssignMags();

public:

  Associations();

  //! Sets a tangent point (reasonably centered for the input image set).
  void SetCommonTangentPoint(const Point &CommonTangentPoint)
    { commonTangentPoint = CommonTangentPoint;};

  //! can be used to project sidereal coordinates related to the image set on a plane.
  Point CommonTangentPoint() const { return commonTangentPoint;}

  //! actually imports the catalog and turns it into a MeasuredStarList using the provided "loader".
//  bool AddImage(const std::string &ReducedImageName);

  //! same as above
//  bool AddImage(const ReducedImage &Ri);
    bool AddImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Ri,
            const PTR(lsst::afw::image::TanWcs) wcs,
            const PTR(lsst::daf::base::PropertySet) meta,
            const lsst::afw::geom::Box2I &bbox,
            const std::string &filter,
            const PTR(lsst::afw::image::Calib) calib  );

  //! incrementaly builds a merged catalog of all image catalogs
  void AssociateCatalogs(const double MatchCutInArcSec = 0, 
			   const bool UseFittedList = false,
			   const bool EnlargeFittedList = true);

  //! refit (independentlyy) all WCSs using the current fittedStarList    
  void RefitWCSs(const unsigned Order);


  //! Collect stars form an external reference catalog (USNO-A by default) that match the FittedStarList. Optionally project these RefStar s on the tangent plane defined by the CommonTangentPoint().
  void CollectRefStars(const bool ProjectOnTP=true);



    

  //! This is Monte-Carlo stuff -- to remove this from associations
  //! the Association class should provide us w/ iterators and acceptors
  //! that way, all the MC stuff could be removed
  void CollectMCStars(int realization = -1);
  //    void CheckMCStars(); // DELETE THIS METHOD
  
  //! This method associates the catalogs with an external 
  //! catalog of photometric ref stars
  void CollectPhotometricRefStars(std::string const& catalogname);
  void AssociatePhotometricRefStars(double MatchCutInArcSec);
  unsigned int GetNbAssociatedPhotometricRefStars() const { return nb_photref_associations; }

  //! Sends back the fitted stars coordinates on the sky FittedStarsList::inTangentPlaneCoordinates keeps track of that.
  void DeprojectFittedStars();

  //! Set the color field of FittedStar 's from a colored catalog. 
  /* If Color is "g-i", then the color is assigned from columns "g" and "i" of the colored catalog. */
  void SetFittedStarColors(std::string DicStarListName, 
				       std::string Color,
			   const double &MatchCutArcSec);
  
  //    void SetRefPhotFactor(int chip, double photfact);
  
  //! apply cuts (mainly number of measurements) on potential FittedStars
  void SelectFittedStars();
  
  void RefineFittedStars();

  void check_refcounts();
  
  void PrintStats() const;
  
  const CcdImageList& TheCcdImageList() const {return ccdImageList;}
  
  void SetFittedStarList(const FittedStarList &L);
  
  void AssignIndexes();
  
  unsigned int NShoots() const { return nshoots_; }

  //! Number of different bands in the input image list. Not implemented so far
  unsigned NBands() const {return 1;}

  
private:
  void AssociateRefStars(const double &MatchCutInArcSec, const Gtransfo *T);
  unsigned int nshoots_;
  unsigned int nb_photref_associations;
};

#endif /* ASSOCIATIONS__H */

}}} // end of namespaces
