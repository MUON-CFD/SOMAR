#ifndef Md2c23a81b8e7d2ba821ba765e1d44345
#error Do not include LevelDataFluxBox.H directly.
#endif

#ifndef USE_LEVELDATAFLUXBOXMOD
#define USE_LEVELDATAFLUXBOXMOD

#include <sstream>
#include "FluxBox.H"
#include "MayDay.H"
#include "StaggeredCopier.H"
#include <array>
#include <memory>

template <>
class LevelData<FluxBox> : public BoxLayoutData<FluxBox>
{
public:
    typedef FluxBox DataHolderT;
    typedef std::array<std::shared_ptr<LDOperator<FArrayBox>>, SpaceDim> LDOpT;
    typedef std::array<StaggeredCopier, SpaceDim> CopierT;

    // -------------------------------------------------------------------------
    /// \{
    /// \name Construction / destruction

    ///
    virtual void
    clear();

    ///
    LevelData();

    ///
    LevelData(
        const DisjointBoxLayout&    dp,
        int                         comps,
        const IntVect&              ghost     = IntVect::Zero,
        const DataFactory<FluxBox>& a_factory = DefaultDataFactory<FluxBox>());

#ifdef CH_USE_THRUST
    /// Special copy constructor: applies only if the argument is a LevelData
    /// defined on Cudified containers
    template <
        class S,
        typename std::enable_if<
            IsCudified<typename S::DataHolderT>::value &&
            std::is_same<S, LevelData<typename S::DataHolderT>>::value>::type* =
            nullptr>
    LevelData(
        const S&                    a_src,
        const DataFactory<FluxBox>& a_factory = DefaultDataFactory<FluxBox>())
    {
        this->define(
            a_src.getBoxes(), a_src.nComp(), a_src.ghostVect(), a_factory);

        auto dit = a_src.getBoxes().dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            auto buf = this->operator[](dit).dataPtr();

            thrust::copy(a_src[dit].begin(), a_src[dit].end(), buf);
        }
    }
#endif

    ///
    virtual ~LevelData();

    ///
    virtual void
    define(
        const DisjointBoxLayout&    dp,
        int                         comps,
        const IntVect&              ghost     = IntVect::Zero,
        const DataFactory<FluxBox>& a_factory = DefaultDataFactory<FluxBox>());

    ///
    /**
      copy definer.  this LevelData thrown away and da's LevelData copied
      */
    virtual void
    define(
        const LevelData<FluxBox>&   da,
        const DataFactory<FluxBox>& a_factory = DefaultDataFactory<FluxBox>());

    ///
    /**
      Copy definer.  'this' LevelData thrown away and da's LevelData copied.
      This LevelData is now defined for the range [0, comps.size()] and filled
      with the data in da from [comps.begin(), comps.end()]
      */
    virtual void
    define(
        const LevelData<FluxBox>&   da,
        const Interval&             comps,
        const DataFactory<FluxBox>& a_factory = DefaultDataFactory<FluxBox>());
    /// \}

    // -------------------------------------------------------------------------
    /// \{
    /// \name LevelData --> LevelData copyTo functions

    /// Simplest case -- assumes source and dest have same interval
    virtual void
    copyTo(LevelData<FluxBox>& a_dest) const;

    /// Same as copyTo that takes a BoxLayoutData, except that it fills the
    /// ghost cells of 'dest' with data from 'this' also.
    virtual void
    copyTo(const Interval&     a_srcComps,
           LevelData<FluxBox>& a_dest,
           const Interval&     a_destComps) const;

    /// Simplest case -- assumes source and dest have same interval
    virtual void
    copyTo(LevelData<FluxBox>& a_dest,
           const CopierT&      a_copier,
           const LDOpT&        a_op = LDOpT()) const;

    /// same as copyTo that takes a BoxLayoutData, except that it fills the
    /// ghost cells of 'dest' with data from 'this' also. USer passes in
    /// a prebuilt Copier object.
    virtual void
    copyTo(const Interval&     a_srcComps,
           LevelData<FluxBox>& a_dest,
           const Interval&     a_destComps,
           const CopierT&      a_copier,
           const LDOpT&        a_op = LDOpT()) const;
    /// \}

    // -------------------------------------------------------------------------
    /// \{
    /// \name LevelData --> BoxLayoutData copyTo functions

    /// Simplest case -- assumes source and dest have same interval
    virtual void
    copyTo(BoxLayoutData<FluxBox>& a_dest) const;

    ///
    virtual void
    copyTo(const Interval&         a_srcComps,
           BoxLayoutData<FluxBox>& a_dest,
           const Interval&         a_destComps) const;

    /// Simplest case -- assumes source and dest have same interval
    virtual void
    copyTo(BoxLayoutData<FluxBox>& a_dest,
           const CopierT&          a_copier,
           const LDOpT&            a_op = LDOpT()) const;

    /// same as copyTo that takes a BoxLayoutData, except that it fills the
    /// ghost cells of 'dest' with data from 'this' also. USer passes in
    /// a prebuilt Copier object
    virtual void
    copyTo(const Interval&         a_srcComps,
           BoxLayoutData<FluxBox>& a_dest,
           const Interval&         a_destComps,
           const CopierT&          a_copier,
           const LDOpT&            a_op = LDOpT()) const;
    /// \}

    // -------------------------------------------------------------------------
    /// \{
    /// \name Exchange functions

    /// Simplest case -- do all components
    virtual void
    exchange(void);

    ///
    virtual void
    exchange(const Interval& comps);

    ///
    virtual void
    exchange(const Interval& a_comps, const CopierT& a_copier);

    /// Simplest case -- do all components
    virtual void
    exchange(const CopierT& a_copier);

    /// asynchronous exchange start.  load and fire off messages.
    virtual void
    exchangeBegin(const CopierT& a_copier);

    /// finish asynchronous exchange
    virtual void
    exchangeEnd();

    ///
    virtual void
    exchangeNoOverlap(const CopierT& a_copier);
    /// \}

    // -------------------------------------------------------------------------
    /// \{
    /// \name Basic accessors

    ///
    inline const IntVect&
    ghostVect() const
    {
        return m_ghost;
    }

    ///
    inline const DisjointBoxLayout&
    getBoxes() const
    {
        return m_disjointBoxLayout;
    }

    ///
    inline const DisjointBoxLayout&
    disjointBoxLayout() const
    {
        return m_disjointBoxLayout;
    }
    /// \}

    // -------------------------------------------------------------------------
    ///
    /** User writes a function with the signature:

    \code
       void myfunction(const Box& box, int n_comps, FluxBox& t)
       {
         your code here;
       }
    \endcode

       They can then hand this off to LayoutData::apply, which invokes this
       function for every FluxBox.  The argument "box" is the Box (as known to
    the DisjointBoxLayout here) associated with that FluxBox and the argument
    "n_comps" is the number of components in this LevelData.

       Your function must not be inline.

    For example:
    \code
    LevelData<FArrayBox> data(layout, 3, IntVect::Unit);
    struct val
    {
      static void set1(const Box& box, int n_comps, const FArrayBox& fab)
      {
        fab.setVal( box.smallEnd(0), box, 0, n_comps );
      }
    };

    data.apply(val::set1);
    \endcode
  */
    virtual void
    apply(void (*a_Function)(const Box&, int, FluxBox&));

    /** For use with apply( const ApplyFunctor& ) */
    struct ApplyFunctor {
        virtual ~ApplyFunctor()
        {
        }

        virtual void
        operator()(const Box&, int, FluxBox&) const = 0;
    };

    /** Like the other apply(), but here the argument is an instance of a class
      derived from LevelData::ApplyFunctor, and which implements ApplyFunctor's
      pure virtual  void operator()(const Box& box, int n_comps, FluxBox& t)
    const.

      Going with an instance of such a class is more convenient if you want
      the thing you pass to apply() to have state.

    For example:
    \code
        class MyFunctor : public LevelData<FArrayBox>::ApplyFunctor
        {
        public:
          MyFunctor( Real x ) : m_x(x)
          {
            ...
          }
          virtual void operator()(const Box& box, int n_comps, FArrayBox& fab)
    const
          {
            fab.setVal( m_x, box, 0, n_comps );
          }
        private:
          const Real m_x;
        }

        LevelData<FArrayBox> data(layout, 3, IntVect::Unit);
        data.apply( MyFunctor(3.14, 0) );

    \endcode

  */
    virtual void
    apply(const ApplyFunctor&);

    void
    degenerate(LevelData<FluxBox>& a_to, const SliceSpec& a_ss) const;

    // -------------------------------------------------------------------------
    /**
       \name overidden virtual functions

      These functions will invoke error messages when invoked.
      C++ will ensure that constructors are not called for the
      base class by a user, but a 'define' function has no such protection,
      hence the need to prevent such usage.  A runtime error is not
      a perfect solution...(strong construction gets around this  *sigh*).
      classes that derive from LevelData will have to turn its valid
      defines into runtime errors also and make its own defines.  Thus
      taking over the job of the compiler.
      */
    /*@{*/
    virtual void
    define(const BoxLayout& dp, int comps, const DataFactory<FluxBox>& factory);

    ///
    virtual void
    define(const BoxLayoutData<FluxBox>& da,
           const DataFactory<FluxBox>& factory = DefaultDataFactory<FluxBox>());

    ///
    virtual void
    define(const BoxLayoutData<FluxBox>& da,
           const Interval&               comps,
           const DataFactory<FluxBox>& factory = DefaultDataFactory<FluxBox>());

    virtual void
    define(const BoxLayout& deadFunction);
    /*@}*/

#ifdef CH_USE_PYTHON
  PyObject* pack() const;
  PyObject* pack();
#endif

protected:
    DisjointBoxLayout m_disjointBoxLayout;

    IntVect m_ghost;

    friend void
    aliasLevelData<FluxBox>(LevelData<FluxBox>& a_alias,
                            LevelData<FluxBox>* a_original,
                            const Interval&     a_interval);

    CopierT m_exCopier;
    CopierT m_exCornerCopier1;
    CopierT m_exCornerCopier2;
    std::array<LevelData<FArrayBox>, CH_SPACEDIM> m_alias;
};


// -----------------------------------------------------------------------------
// These functions alias an interval of a LevelData<FluxBox>.
// -----------------------------------------------------------------------------

// This is modeled after aliasLevelData, which takes a pointer for the
// second argument.
#if defined(CH_USE_THRUST) || defined(__NVCC__)
    template <typename M> void
    aliasFluxBox(CudaFluxBox<M>& a_new, CudaFluxBox<M>* a_oldPtr, const Interval& a_ivl);
#else
    void
    aliasFluxBox(FluxBox& a_new, FluxBox* a_oldPtr, const Interval& a_ivl);
#endif

// This can be passed into the LevelData constructor to create an alias.
template <>
class AliasDataFactory<FluxBox> : public DataFactory<FluxBox>
{
public:
    AliasDataFactory(BoxLayoutData<FluxBox>* a_origPtr, const Interval& a_ivl)
    : m_origPtr(a_origPtr)
    , m_ivl(a_ivl)
    {
        // TODO: Sanity checks.
    }

    virtual ~AliasDataFactory()
    {}

    // factory function.  creates a new 'T' object using an aliased dataPtr for
    // T creates a new 'T' object and returns a pointer to it.  Responsiblitly
    // for calling operator 'delete' on this pointer is passed to the user.
    virtual FluxBox*
    create(const Box&       /*a_box*/,
           int              a_nComps [[maybe_unused]],
           const DataIndex& a_di) const
    {
        CH_assert(a_nComps = m_ivl.size());
        FluxBox* newPtr = new FluxBox;
        aliasFluxBox(*newPtr, &((*m_origPtr)[a_di]), m_ivl);
        return newPtr;
    }

protected:
    BoxLayoutData<FluxBox>* m_origPtr;
    Interval                m_ivl;
};


// The LevelData<FluxBox> version of aliasLevelData.
template <>
void
aliasLevelData<FluxBox>(LevelData<FluxBox>& a_alias,
                        LevelData<FluxBox>* a_original,
                        const Interval&     a_interval);


// -----------------------------------------------------------------------------
// These functions alias an the normal components LevelData<FluxBox>.
// -----------------------------------------------------------------------------

// This is modeled after aliasLevelData, which takes a pointer for the
// second argument. *a_oldPtr must have SpaceDim comps.
#if defined(CH_USE_THRUST) || defined(__NVCC__)
    template <typename M> void
    aliasNormalFluxes(CudaFluxBox<M>& a_new, CudaFluxBox<M>* a_oldPtr);
#else
    void aliasNormalFluxes(FluxBox& a_new, FluxBox* a_oldPtr);
#endif

// This can be passed into the LevelData constructor to create an alias.
class AliasNormalDataFactory: public DataFactory<FluxBox>
{
public:
    AliasNormalDataFactory(BoxLayoutData<FluxBox>* a_origPtr);

    virtual ~AliasNormalDataFactory();

    // factory function.  creates a new 'T' object using an aliased dataPtr for
    // T creates a new 'T' object and returns a pointer to it.  Responsiblitly
    // for calling operator 'delete' on this pointer is passed to the user.
    virtual FluxBox*
    create(const Box& a_box, int a_nComps, const DataIndex& a_di) const;

protected:
    BoxLayoutData<FluxBox>* m_origPtr;
};


// Similar to aliasLevelData. *a_original must have SpaceDim comps.
#if defined(CH_USE_THRUST) || defined(__NVCC__)
    template <typename M> void
    aliasNormalFluxes(LevelData<CudaFluxBox<M>>& a_alias,
                    LevelData<CudaFluxBox<M>>* a_original);
#else
    void
    aliasNormalFluxes(LevelData<FluxBox>& a_alias,
                    LevelData<FluxBox>* a_original);
#endif


#endif  // USE_LEVELDATAFLUXBOXMOD
