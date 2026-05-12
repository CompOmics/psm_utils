"""
SQLAlchemy ORM models for Proteome Discoverer MSF database files.

This module provides SQLAlchemy table definitions for interfacing with Proteome Discoverer MSF
(Mascot Search Form) database files. MSF files contain proteomics search results including
peptide identifications, protein annotations, spectra metadata, and quantification data.

The table definitions are auto-generated from MSF schema and follow SQLAlchemy 2.0 patterns
with proper typing support.

Examples
--------
>>> from psm_utils.io._pd_msf_tables import Base, Peptide
>>> # Use with SQLAlchemy session to query MSF database
>>> session.query(Peptide).filter(Peptide.ConfidenceLevel > 2).all()

Notes
-----
These models are primarily used internally by the proteome_discoverer module for reading PSM
data from MSF files.

"""

from __future__ import annotations

from datetime import datetime

from sqlalchemy import (
    CHAR,
    BigInteger,
    Boolean,
    DateTime,
    Float,
    Index,
    Integer,
    LargeBinary,
    MetaData,
    SmallInteger,
    String,
    Text,
    UniqueConstraint,
    text,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column


class Base(DeclarativeBase):
    """Base class for all MSF table models."""

    pass


# Module-level metadata reference for table definitions
metadata: MetaData = Base.metadata


class AminoAcidModification(Base):
    __tablename__ = "AminoAcidModifications"

    AminoAcidModificationID: Mapped[int] = mapped_column(Integer, primary_key=True)
    ModificationName: Mapped[str] = mapped_column(String, nullable=False)
    DeltaMass: Mapped[float | None] = mapped_column(Float)
    Substitution: Mapped[str | None] = mapped_column(String)
    LeavingGroup: Mapped[str | None] = mapped_column(String)
    Abbreviation: Mapped[str] = mapped_column(String, nullable=False)
    PositionType: Mapped[int] = mapped_column(Integer, nullable=False)
    IsActive: Mapped[bool | None] = mapped_column(Boolean)
    DeltaAverageMass: Mapped[float | None] = mapped_column(Float)
    UnimodAccession: Mapped[str | None] = mapped_column(String)
    IsSubstitution: Mapped[bool] = mapped_column(Boolean, nullable=False, server_default=text("0"))


class AminoAcidModificationsAminoAcid(Base):
    __tablename__ = "AminoAcidModificationsAminoAcids"

    AminoAcidModificationID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    AminoAcidID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Classification: Mapped[int] = mapped_column(Integer, nullable=False)


class AminoAcidModificationsAminoAcidsNL(Base):
    __tablename__ = "AminoAcidModificationsAminoAcidsNL"

    AminoAcidModificationID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    AminoAcidID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    NeutralLossID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)


class AminoAcidModificationsNeutralLoss(Base):
    __tablename__ = "AminoAcidModificationsNeutralLosses"

    NeutralLossID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Name: Mapped[str] = mapped_column(String, nullable=False)
    MonoisotopicMass: Mapped[float] = mapped_column(Float, nullable=False)
    AverageMass: Mapped[float] = mapped_column(Float, nullable=False)


class AminoAcid(Base):
    __tablename__ = "AminoAcids"

    AminoAcidID: Mapped[int] = mapped_column(Integer, primary_key=True)
    AminoAcidName: Mapped[str] = mapped_column(String, nullable=False)
    OneLetterCode: Mapped[str | None] = mapped_column(CHAR)
    ThreeLetterCode: Mapped[str | None] = mapped_column(CHAR)
    MonoisotopicMass: Mapped[float] = mapped_column(Float, nullable=False)
    AverageMass: Mapped[float] = mapped_column(Float, nullable=False)
    SumFormula: Mapped[str | None] = mapped_column(String)


class AnnotationDataVersion(Base):
    __tablename__ = "AnnotationDataVersion"

    PcDataVersion: Mapped[int] = mapped_column(Integer, primary_key=True)
    PcDataRelease: Mapped[int] = mapped_column(BigInteger, nullable=False)


class AnnotationDataset(Base):
    __tablename__ = "AnnotationDataset"

    DatasetId: Mapped[int] = mapped_column(Integer, primary_key=True)
    Name: Mapped[str] = mapped_column(String, nullable=False)
    DisplayName: Mapped[str] = mapped_column(String, nullable=False)
    Guid: Mapped[str] = mapped_column(String, nullable=False)
    Description: Mapped[str | None] = mapped_column(Text)


class AnnotationGroup(Base):
    __tablename__ = "AnnotationGroups"

    AnnotationGroupId: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Description: Mapped[str | None] = mapped_column(Text)
    DatasetId: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Position: Mapped[int] = mapped_column(Integer, nullable=False)
    ColorR: Mapped[int] = mapped_column(Integer, nullable=False)
    ColorG: Mapped[int] = mapped_column(Integer, nullable=False)
    ColorB: Mapped[int] = mapped_column(Integer, nullable=False)
    GroupDefinition: Mapped[bytes | None] = mapped_column(LargeBinary)


class AnnotationType(Base):
    __tablename__ = "AnnotationTypes"

    AnnotationTypeId: Mapped[int] = mapped_column(Integer, primary_key=True)
    Name: Mapped[str] = mapped_column(String, nullable=False)
    Description: Mapped[str | None] = mapped_column(Text)


class Annotation(Base):
    __tablename__ = "Annotations"

    AnnotationId: Mapped[int] = mapped_column(Integer, primary_key=True)
    Accession: Mapped[str] = mapped_column(String, nullable=False)
    Description: Mapped[str | None] = mapped_column(Text)
    type: Mapped[int | None] = mapped_column(Integer)


class AnnotationsAnnotationGroup(Base):
    __tablename__ = "AnnotationsAnnotationGroups"

    AnnotationId: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    AnnotationGroupId: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)


class AnnotationsProtein(Base):
    __tablename__ = "AnnotationsProtein"

    proteinID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    AnnotationId: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Evidence: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    PositionBegin: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    PositionEnd: Mapped[int | None] = mapped_column(Integer)
    ProteinAccession: Mapped[str] = mapped_column(String, primary_key=True, nullable=False)


class Chromatogram(Base):
    __tablename__ = "Chromatograms"

    FileID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    TraceType: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Chromatogram: Mapped[str] = mapped_column(String, nullable=False)


class CustomDataField(Base):
    __tablename__ = "CustomDataFields"

    FieldID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Guid: Mapped[str] = mapped_column(String, nullable=False)
    DisplayName: Mapped[str] = mapped_column(String, nullable=False)
    SourceNodeNumber: Mapped[int] = mapped_column(Integer, nullable=False)
    TargetNodeNumber: Mapped[int] = mapped_column(Integer, nullable=False)
    DataType: Mapped[int] = mapped_column(Integer, nullable=False)
    DataTarget: Mapped[int] = mapped_column(Integer, nullable=False)
    Version: Mapped[float] = mapped_column(Float, nullable=False)
    AccessMode: Mapped[int | None] = mapped_column(Integer, server_default=text("0"))
    Visibility: Mapped[int | None] = mapped_column(Integer, server_default=text("0"))
    GroupVisibility: Mapped[int | None] = mapped_column(Integer, server_default=text("0"))
    Format: Mapped[str | None] = mapped_column(String)
    PlotType: Mapped[int] = mapped_column(Integer, nullable=False)
    DataPurpose: Mapped[str | None] = mapped_column(String)


class CustomDataPeptide(Base):
    __tablename__ = "CustomDataPeptides"

    FieldID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue: Mapped[str | None] = mapped_column(String)


class CustomDataPeptidesDecoy(Base):
    __tablename__ = "CustomDataPeptides_decoy"

    FieldID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue: Mapped[str | None] = mapped_column(String)


class CustomDataProcessingNode(Base):
    __tablename__ = "CustomDataProcessingNodes"

    FieldID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ProcessingNodeNumber: Mapped[int] = mapped_column(
        Integer, primary_key=True, nullable=False, index=True
    )
    FieldValue: Mapped[str | None] = mapped_column(String)


class CustomDataProtein(Base):
    __tablename__ = "CustomDataProteins"

    FieldID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ProteinID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue: Mapped[str | None] = mapped_column(String)


class CustomDataProteinsDecoy(Base):
    __tablename__ = "CustomDataProteins_decoy"

    FieldID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ProteinID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue: Mapped[str | None] = mapped_column(String)


class CustomDataSpectra(Base):
    __tablename__ = "CustomDataSpectra"

    FieldID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    SpectrumID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue: Mapped[str | None] = mapped_column(String)


class Enzyme(Base):
    __tablename__ = "Enzymes"

    EnzymeID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Name: Mapped[str] = mapped_column(String, nullable=False)
    Abbreviation: Mapped[str] = mapped_column(String, nullable=False)
    Seperator: Mapped[str] = mapped_column(String, nullable=False)
    NonSeperator: Mapped[str] = mapped_column(String, nullable=False)
    Offset: Mapped[int] = mapped_column(Integer, nullable=False)


class EnzymesCleavageSpecificity(Base):
    __tablename__ = "EnzymesCleavageSpecificities"

    EnzymeID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Specificity: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)


class EventAnnotation(Base):
    __tablename__ = "EventAnnotations"
    __table_args__ = (
        Index(
            "IX_EventAnnotations_IsotopePatternID_QuanResultID", "IsotopePatternID", "QuanResultID"
        ),
        Index("IX_EventAnnotations_QuanResultID_QuanChannelID", "QuanResultID", "QuanChannelID"),
    )

    EventID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Charge: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    IsotopePatternID: Mapped[int] = mapped_column(Integer, nullable=False)
    QuanResultID: Mapped[int] = mapped_column(Integer, nullable=False)
    QuanChannelID: Mapped[int] = mapped_column(Integer, nullable=False)


class EventAreaAnnotation(Base):
    __tablename__ = "EventAreaAnnotations"

    EventID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Charge: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    IsotopePatternID: Mapped[int] = mapped_column(Integer, nullable=False, index=True)
    QuanResultID: Mapped[int] = mapped_column(Integer, nullable=False)


class Event(Base):
    __tablename__ = "Events"
    __table_args__ = (
        Index("IX_Events_FileID_LeftRT_RightRT", "FileID", "LeftRT", "RightRT"),
        Index("IX_Events_FileID_RT", "FileID", "RT"),
    )

    EventID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Mass: Mapped[float] = mapped_column(Float, nullable=False)
    MassAvg: Mapped[float] = mapped_column(Float, nullable=False)
    Area: Mapped[float] = mapped_column(Float, nullable=False)
    Intensity: Mapped[float] = mapped_column(Float, nullable=False)
    PeakWidth: Mapped[float] = mapped_column(Float, nullable=False)
    RT: Mapped[float] = mapped_column(Float, nullable=False)
    LeftRT: Mapped[float] = mapped_column(Float, nullable=False)
    RightRT: Mapped[float] = mapped_column(Float, nullable=False)
    SN: Mapped[float] = mapped_column(Float, nullable=False, server_default=text("0.0"))
    FileID: Mapped[int] = mapped_column(Integer, nullable=False)


class FastaFile(Base):
    __tablename__ = "FastaFiles"

    FastaFileID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    FileName: Mapped[str] = mapped_column(String, nullable=False)
    State: Mapped[int] = mapped_column(Integer, nullable=False)
    VirtualFileName: Mapped[str] = mapped_column(String, nullable=False)
    FileSize: Mapped[int] = mapped_column(BigInteger, nullable=False)
    FileTime: Mapped[int] = mapped_column(BigInteger, nullable=False)
    NumberOfProteins: Mapped[int | None] = mapped_column(BigInteger)
    NumberOfAminoAcids: Mapped[int | None] = mapped_column(BigInteger)
    FileHashCode: Mapped[int | None] = mapped_column(BigInteger)
    Hidden: Mapped[bool] = mapped_column(Boolean, nullable=False)
    IsSrfImport: Mapped[bool] = mapped_column(Boolean, nullable=False)
    IsScheduledForDeletion: Mapped[bool] = mapped_column(
        Boolean, nullable=False, server_default=text("0")
    )


class FastaFilesProteinAnnotation(Base):
    __tablename__ = "FastaFilesProteinAnnotations"

    FastaFileID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ProteinAnnotationID: Mapped[int] = mapped_column(
        Integer, primary_key=True, nullable=False, index=True
    )


class FileInfo(Base):
    __tablename__ = "FileInfos"

    FileID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    FileName: Mapped[str] = mapped_column(String, nullable=False)
    FileTime: Mapped[str] = mapped_column(String, nullable=False)
    FileSize: Mapped[int] = mapped_column(BigInteger, nullable=False)
    PhysicalFileName: Mapped[str] = mapped_column(String, nullable=False)
    FileType: Mapped[int] = mapped_column(SmallInteger, nullable=False)


class MassPeakRelation(Base):
    __tablename__ = "MassPeakRelations"

    MassPeakID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    RelatedMassPeakID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)


class MassPeak(Base):
    __tablename__ = "MassPeaks"

    MassPeakID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    Charge: Mapped[int | None] = mapped_column(SmallInteger)
    Intensity: Mapped[float | None] = mapped_column(Float)
    Mass: Mapped[float | None] = mapped_column(Float)
    ScanNumbers: Mapped[str | None] = mapped_column(String)
    FileID: Mapped[int | None] = mapped_column(Integer)
    PercentIsolationInterference: Mapped[float | None] = mapped_column(Float)
    IonInjectTime: Mapped[int | None] = mapped_column(Integer)


class PeptideScore(Base):
    __tablename__ = "PeptideScores"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ScoreID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ProcessingNodeID: Mapped[int | None] = mapped_column(Integer)
    ScoreValue: Mapped[float] = mapped_column(Float, nullable=False)


class PeptideScoreDecoy(Base):
    __tablename__ = "PeptideScores_decoy"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ScoreID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ProcessingNodeID: Mapped[int | None] = mapped_column(Integer)
    ScoreValue: Mapped[float] = mapped_column(Float, nullable=False)


class Peptide(Base):
    __tablename__ = "Peptides"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False, index=True)
    SpectrumID: Mapped[int] = mapped_column(Integer, nullable=False, index=True)
    TotalIonsCount: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    MatchedIonsCount: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    ConfidenceLevel: Mapped[int] = mapped_column(SmallInteger, nullable=False)
    SearchEngineRank: Mapped[int] = mapped_column(Integer, nullable=False)
    Hidden: Mapped[bool] = mapped_column(Boolean, nullable=False, server_default=text("0"))
    Sequence: Mapped[str | None] = mapped_column(String)
    Annotation: Mapped[str | None] = mapped_column(String)
    UniquePeptideSequenceID: Mapped[int] = mapped_column(
        Integer, nullable=False, server_default=text("1")
    )
    MissedCleavages: Mapped[int] = mapped_column(SmallInteger, nullable=False)


class PeptidesAminoAcidModification(Base):
    __tablename__ = "PeptidesAminoAcidModifications"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True)
    AminoAcidModificationID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Position: Mapped[int] = mapped_column(Integer, primary_key=True)


class PeptidesAminoAcidModificationsDecoy(Base):
    __tablename__ = "PeptidesAminoAcidModifications_decoy"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True)
    AminoAcidModificationID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Position: Mapped[int] = mapped_column(Integer, primary_key=True)


class PeptidesProtein(Base):
    __tablename__ = "PeptidesProteins"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True, index=True)
    ProteinID: Mapped[int] = mapped_column(Integer, primary_key=True)


class PeptidesProteinDecoy(Base):
    __tablename__ = "PeptidesProteins_decoy"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True, index=True)
    ProteinID: Mapped[int] = mapped_column(Integer, primary_key=True)


class PeptidesReferenceSpectra(Base):
    __tablename__ = "PeptidesReferenceSpectra"

    PeptideID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    ReferenceSpectrumID: Mapped[int] = mapped_column(Integer)


class PeptidesTerminalModification(Base):
    __tablename__ = "PeptidesTerminalModifications"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True)
    TerminalModificationID: Mapped[int] = mapped_column(Integer, primary_key=True)


class PeptidesTerminalModificationDecoy(Base):
    __tablename__ = "PeptidesTerminalModifications_decoy"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True)
    TerminalModificationID: Mapped[int] = mapped_column(Integer, primary_key=True)


class PeptideDecoy(Base):
    __tablename__ = "Peptides_decoy"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    PeptideID: Mapped[int] = mapped_column(Integer, primary_key=True, index=True)
    SpectrumID: Mapped[int] = mapped_column(Integer, index=True)
    TotalIonsCount: Mapped[int] = mapped_column(SmallInteger)
    MatchedIonsCount: Mapped[int] = mapped_column(SmallInteger)
    ConfidenceLevel: Mapped[int] = mapped_column(SmallInteger)
    SearchEngineRank: Mapped[int] = mapped_column(Integer)
    Sequence: Mapped[str | None] = mapped_column(String)
    Annotation: Mapped[str | None] = mapped_column(String)
    UniquePeptideSequenceID: Mapped[int] = mapped_column(Integer, server_default=text("1"))
    MissedCleavages: Mapped[int] = mapped_column(SmallInteger)


class PrecursorIonAreaSearchSpectra(Base):
    __tablename__ = "PrecursorIonAreaSearchSpectra"

    QuanResultID: Mapped[int] = mapped_column(
        Integer, primary_key=True, nullable=False, index=True
    )
    SearchSpectrumID: Mapped[int | None] = mapped_column(Integer, primary_key=True)


class PrecursorIonQuanResult(Base):
    __tablename__ = "PrecursorIonQuanResults"
    __table_args__ = (
        Index(
            "IX_PrecursorIonQuanResults_QuanResultID_QuanChannelID",
            "QuanResultID",
            "QuanChannelID",
        ),
    )

    QuanChannelID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    QuanResultID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Mass: Mapped[float] = mapped_column(Float, nullable=False)
    Charge: Mapped[int] = mapped_column(Integer, nullable=False)
    Area: Mapped[float | None] = mapped_column(Float)
    RetentionTime: Mapped[float | None] = mapped_column(Float)


class PrecursorIonQuanResultsSearchSpectra(Base):
    __tablename__ = "PrecursorIonQuanResultsSearchSpectra"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    QuanResultID: Mapped[int] = mapped_column(
        Integer, primary_key=True, nullable=False, index=True
    )
    SearchSpectrumID: Mapped[int | None] = mapped_column(Integer, index=True)


class ProcessingNodeConnectionPoint(Base):
    __tablename__ = "ProcessingNodeConnectionPoints"

    ProcessingNodeID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Interface: Mapped[str] = mapped_column(String, primary_key=True, nullable=False)
    ConnectionDirection: Mapped[int] = mapped_column(Integer, nullable=False)
    ConnectionMode: Mapped[int] = mapped_column(Integer, nullable=False)
    ConnectionMultiplicity: Mapped[int] = mapped_column(Integer, nullable=False)
    ConnectionRequirement: Mapped[int] = mapped_column(Integer, nullable=False)
    DataTypeSpecialization: Mapped[str] = mapped_column(String, nullable=False)
    ConnectionDisplayName: Mapped[str] = mapped_column(String, nullable=False)


class ProcessingNodeExtension(Base):
    __tablename__ = "ProcessingNodeExtensions"

    ExtensionID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer)
    Guid: Mapped[str] = mapped_column(String)
    Purpose: Mapped[str] = mapped_column(String)
    PurposeDetail: Mapped[str | None] = mapped_column(String)
    MajorVersion: Mapped[int] = mapped_column(Integer)
    MinorVersion: Mapped[int] = mapped_column(Integer)
    Settings: Mapped[str | None] = mapped_column(Text)


class ProcessingNodeFilterParameter(Base):
    __tablename__ = "ProcessingNodeFilterParameters"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    FilterParameterName: Mapped[str] = mapped_column(String, primary_key=True)
    FilterModuleTypeID: Mapped[int] = mapped_column(Integer)
    FilterModuleNumber: Mapped[int] = mapped_column(Integer)
    ProcessingNodeID: Mapped[int] = mapped_column(Integer)
    FilterParameterValue: Mapped[float] = mapped_column(Float)


class ProcessingNodeInterface(Base):
    __tablename__ = "ProcessingNodeInterfaces"

    ProcessingNodeID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    InterfaceKind: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    InterfaceName: Mapped[str] = mapped_column(String, primary_key=True, nullable=False)


class ProcessingNodeParameter(Base):
    __tablename__ = "ProcessingNodeParameters"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    ParameterName: Mapped[str] = mapped_column(String, primary_key=True)
    FriendlyName: Mapped[str] = mapped_column(String)
    ProcessingNodeID: Mapped[int] = mapped_column(Integer)
    IntendedPurpose: Mapped[int] = mapped_column(Integer)
    PurposeDetails: Mapped[str] = mapped_column(String)
    Hidden: Mapped[bool] = mapped_column(Boolean)
    Advanced: Mapped[bool] = mapped_column(Boolean)
    Category: Mapped[str] = mapped_column(String)
    Position: Mapped[int] = mapped_column(Integer)
    ParameterValue: Mapped[str] = mapped_column(String)
    ValueDisplayString: Mapped[str] = mapped_column(String)


class ProcessingNodeScore(Base):
    __tablename__ = "ProcessingNodeScores"
    __table_args__ = (UniqueConstraint("ProcessingNodeID", "ScoreName"),)

    ProcessingNodeID: Mapped[int] = mapped_column(Integer)
    ScoreID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    ScoreName: Mapped[str] = mapped_column(String)
    FriendlyName: Mapped[str] = mapped_column(String)
    Description: Mapped[str] = mapped_column(String)
    FormatString: Mapped[str] = mapped_column(String)
    ScoreCategory: Mapped[int] = mapped_column(Integer)
    Hidden: Mapped[bool] = mapped_column(Boolean)
    IsMainScore: Mapped[bool] = mapped_column(Boolean)
    ScoreGUID: Mapped[str] = mapped_column(String)


class ProcessingNode(Base):
    __tablename__ = "ProcessingNodes"

    ProcessingNodeNumber: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    ProcessingNodeID: Mapped[int] = mapped_column(Integer)
    ProcessingNodeParentNumber: Mapped[str] = mapped_column(String)
    NodeName: Mapped[str | None] = mapped_column(String)
    FriendlyName: Mapped[str] = mapped_column(String)
    MajorVersion: Mapped[int] = mapped_column(Integer)
    MinorVersion: Mapped[int] = mapped_column(Integer)
    NodeComment: Mapped[str | None] = mapped_column(String)
    NodeGUID: Mapped[str] = mapped_column(String)
    ProcessingNodeState: Mapped[int] = mapped_column(Integer, server_default=text("0"))


class ProcessingNodesSpectra(Base):
    __tablename__ = "ProcessingNodesSpectra"

    SendingProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    SpectrumID: Mapped[int] = mapped_column(Integer, primary_key=True, index=True)


class ProteinAnnotation(Base):
    __tablename__ = "ProteinAnnotations"
    __table_args__ = (
        Index(
            "IX_ProteinAnnotations_ProteinID_DescriptionHashCode",
            "ProteinID",
            "DescriptionHashCode",
        ),
    )

    ProteinAnnotationID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    ProteinID: Mapped[int] = mapped_column(Integer)
    DescriptionHashCode: Mapped[int] = mapped_column(BigInteger)
    Description: Mapped[str] = mapped_column(Text)
    TaxonomyID: Mapped[int] = mapped_column(Integer, index=True)


class ProteinIdentificationGroup(Base):
    __tablename__ = "ProteinIdentificationGroups"

    ProteinIdentificationGroupId: Mapped[int] = mapped_column(Integer, primary_key=True)
    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)


class ProteinScore(Base):
    __tablename__ = "ProteinScores"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    ProteinID: Mapped[int] = mapped_column(Integer, primary_key=True)
    ProteinIdentificationGroupID: Mapped[int] = mapped_column(Integer)
    ProteinScore: Mapped[float] = mapped_column(Float)
    Coverage: Mapped[float] = mapped_column(Float, server_default=text("0"))


class ProteinScoresDecoy(Base):
    __tablename__ = "ProteinScores_decoy"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    ProteinID: Mapped[int] = mapped_column(Integer, primary_key=True)
    ProteinIdentificationGroupID: Mapped[int] = mapped_column(Integer)
    ProteinScore: Mapped[float] = mapped_column(Float)
    Coverage: Mapped[float] = mapped_column(Float, server_default=text("0"))


class Protein(Base):
    __tablename__ = "Proteins"

    ProteinID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    Sequence: Mapped[str] = mapped_column(Text)
    SequenceHashCode: Mapped[int] = mapped_column(BigInteger, index=True)
    IsMasterProtein: Mapped[bool] = mapped_column(Boolean, server_default=text("0"))


class ProteinsProteinGroup(Base):
    __tablename__ = "ProteinsProteinGroups"

    ProteinID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    ProteinGroupID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)


class PtmAnnotationDatum(Base):
    __tablename__ = "PtmAnnotationData"

    AnnotationType: Mapped[int] = mapped_column(Integer, primary_key=True)
    ProteinId: Mapped[int] = mapped_column(Integer, primary_key=True)
    AnnotationId: Mapped[int] = mapped_column(Integer, primary_key=True)
    Position: Mapped[int] = mapped_column(Integer, primary_key=True)
    Annotation: Mapped[str | None] = mapped_column(String)


class ReferenceSpectra(Base):
    __tablename__ = "ReferenceSpectra"

    ReferenceSpectrumId: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    Sequence: Mapped[str] = mapped_column(String)
    SequenceHashCode: Mapped[int] = mapped_column(BigInteger)
    Spectrum: Mapped[str] = mapped_column(String)
    SpectrumHashCode: Mapped[int] = mapped_column(BigInteger)
    Comment: Mapped[str | None] = mapped_column(Text)
    CommentHashCode: Mapped[int] = mapped_column(BigInteger)


class ReporterIonQuanResult(Base):
    __tablename__ = "ReporterIonQuanResults"
    __table_args__ = (
        Index(
            "IX_ReporterIonQuanResults_ProcessingNodeNumber_SpectrumID",
            "ProcessingNodeNumber",
            "SpectrumID",
        ),
    )

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    QuanChannelID: Mapped[int] = mapped_column(Integer, primary_key=True)
    SpectrumID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Mass: Mapped[float] = mapped_column(Float)
    Height: Mapped[float | None] = mapped_column(Float)


class ReporterIonQuanResultsSearchSpectra(Base):
    __tablename__ = "ReporterIonQuanResultsSearchSpectra"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    SpectrumID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    SearchSpectrumID: Mapped[int | None] = mapped_column(Integer, index=True)


class ScanEvent(Base):
    __tablename__ = "ScanEvents"

    ScanEventID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    MSLevel: Mapped[int] = mapped_column(Integer)
    Polarity: Mapped[int] = mapped_column(Integer)
    ScanType: Mapped[int] = mapped_column(Integer)
    Ionization: Mapped[int] = mapped_column(Integer)
    MassAnalyzer: Mapped[int] = mapped_column(Integer)
    ActivationType: Mapped[int] = mapped_column(Integer)


class SchemaInfo(Base):
    __tablename__ = "SchemaInfo"

    Version: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    Kind: Mapped[str] = mapped_column(String)
    Date: Mapped[datetime] = mapped_column(DateTime)
    SoftwareVersion: Mapped[str] = mapped_column(String)
    Comment: Mapped[str] = mapped_column(Text)


class Spectrum(Base):
    __tablename__ = "Spectra"

    UniqueSpectrumID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    Spectrum: Mapped[str] = mapped_column(String)
    SpectrumHashCode: Mapped[int | None] = mapped_column(BigInteger)


class SpectrumHeader(Base):
    __tablename__ = "SpectrumHeaders"

    SpectrumID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    MassPeakID: Mapped[int | None] = mapped_column(Integer)
    ScanEventID: Mapped[int | None] = mapped_column(Integer)
    LastScan: Mapped[int | None] = mapped_column(Integer)
    FirstScan: Mapped[int | None] = mapped_column(Integer)
    RetentionTime: Mapped[float | None] = mapped_column(Float)
    Hidden: Mapped[bool] = mapped_column(Boolean, server_default=text("0"))
    ScanNumbers: Mapped[str | None] = mapped_column(String)
    Charge: Mapped[int | None] = mapped_column(SmallInteger)
    Mass: Mapped[float | None] = mapped_column(Float)
    CreatingProcessingNodeNumber: Mapped[int] = mapped_column(Integer)
    UniqueSpectrumID: Mapped[int] = mapped_column(Integer, server_default=text("0"))


class SpectrumScore(Base):
    __tablename__ = "SpectrumScores"

    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer, primary_key=True)
    SpectrumID: Mapped[int] = mapped_column(Integer, primary_key=True)
    Score: Mapped[float] = mapped_column(Float)


class TaxonomyName(Base):
    __tablename__ = "TaxonomyNames"

    TaxonomyID: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False, index=True)
    NameCategory: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    Name: Mapped[str | None] = mapped_column(String)


class TaxonomyNode(Base):
    __tablename__ = "TaxonomyNodes"
    __table_args__ = (
        Index("IX_TaxonomyNodes_LeftNodeIndex_RightNodeIndex", "LeftNodeIndex", "RightNodeIndex"),
    )

    TaxonomyID: Mapped[int | None] = mapped_column(Integer, primary_key=True, unique=True)
    ParentTaxonomyID: Mapped[int] = mapped_column(Integer)
    TaxonomyRank: Mapped[int] = mapped_column(Integer)
    LeftNodeIndex: Mapped[int] = mapped_column(Integer)
    RightNodeIndex: Mapped[int] = mapped_column(Integer)


# TODO: Check which is primary key
class WorkflowInfo(Base):
    __tablename__ = "WorkflowInfo"

    WorkflowGUID: Mapped[str] = mapped_column(String, primary_key=True, nullable=False)
    WorkflowName: Mapped[str] = mapped_column(String, nullable=False)
    WorkflowDescription: Mapped[str] = mapped_column(String, nullable=False)
    WorkflowState: Mapped[int] = mapped_column(Integer, nullable=False, server_default=text("0"))
    WorkflowStartDate: Mapped[datetime] = mapped_column(DateTime, nullable=False)
    WorkflowTemplate: Mapped[str] = mapped_column(String, nullable=False)
    User: Mapped[str] = mapped_column(String, nullable=False)
    MachineGUID: Mapped[str] = mapped_column(String, nullable=False)
    MachineName: Mapped[str] = mapped_column(String, nullable=False)
    MergeSimilarIdentificationResults: Mapped[bool] = mapped_column(Boolean, nullable=False)
    IsValid: Mapped[bool] = mapped_column(Boolean, nullable=False)
    Version: Mapped[int] = mapped_column(Integer, nullable=False)


class WorkflowMessage(Base):
    __tablename__ = "WorkflowMessages"

    MessageID: Mapped[int | None] = mapped_column(Integer, primary_key=True)
    ProcessingNodeID: Mapped[int] = mapped_column(Integer)
    ProcessingNodeNumber: Mapped[int] = mapped_column(Integer)
    Time: Mapped[int] = mapped_column(BigInteger)
    MessageKind: Mapped[int] = mapped_column(Integer)
    Message: Mapped[str] = mapped_column(String)
